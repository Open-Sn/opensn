// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/volume_postprocessor.h"
#include "framework/parameters/input_parameters.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>

namespace opensn
{

namespace
{

/**
 * Uniform grid of bins over the centroids of the local cells. It lists the cells whose centroid
 * can lie in an axis-aligned box without testing every local cell.
 */
class CentroidBins
{
public:
  explicit CentroidBins(const Mesh& grid)
  {
    const auto& cells = grid.GetLocalCells();
    const auto num_cells = cells.size();
    if (num_cells == 0)
      return;

    std::array<double, 3> cmax{};
    for (std::size_t d = 0; d < 3; ++d)
    {
      origin_[d] = std::numeric_limits<double>::infinity();
      cmax[d] = -std::numeric_limits<double>::infinity();
    }
    for (const auto& cell : cells)
      for (std::size_t d = 0; d < 3; ++d)
      {
        origin_[d] = std::min(origin_[d], cell->centroid[d]);
        cmax[d] = std::max(cmax[d], cell->centroid[d]);
      }

    // About one bin per cell, with bins of equal width in the directions the centroids span. A
    // direction thinner than one bin gets a single bin, and the width is recomputed without it.
    std::array<bool, 3> active{};
    for (std::size_t d = 0; d < 3; ++d)
      active[d] = cmax[d] > origin_[d];
    double width = 0.0;
    for (bool changed = true; changed;)
    {
      changed = false;
      double measure = 1.0;
      int num_active = 0;
      for (std::size_t d = 0; d < 3; ++d)
        if (active[d])
        {
          measure *= cmax[d] - origin_[d];
          ++num_active;
        }
      if (num_active == 0)
        break;
      width = std::pow(measure / static_cast<double>(num_cells), 1.0 / num_active);
      for (std::size_t d = 0; d < 3; ++d)
        if (active[d] and cmax[d] - origin_[d] < width)
        {
          active[d] = false;
          changed = true;
        }
    }
    for (std::size_t d = 0; d < 3; ++d)
    {
      num_bins_[d] = 1;
      if (active[d])
        num_bins_[d] = static_cast<std::size_t>(std::ceil((cmax[d] - origin_[d]) / width));
      num_bins_[d] = std::clamp<std::size_t>(num_bins_[d], 1, num_cells);
      bin_width_[d] = (cmax[d] - origin_[d]) / static_cast<double>(num_bins_[d]);
    }

    // Cell IDs grouped by bin, in increasing order within each bin.
    const auto num_bins = num_bins_[0] * num_bins_[1] * num_bins_[2];
    std::vector<std::size_t> cell_bin(num_cells);
    offsets_.assign(num_bins + 1, 0);
    for (const auto& cell : cells)
    {
      const auto& c = cell->centroid;
      cell_bin[cell->local_id] =
        (BinIndex(0, c.x) * num_bins_[1] + BinIndex(1, c.y)) * num_bins_[2] + BinIndex(2, c.z);
      ++offsets_[cell_bin[cell->local_id] + 1];
    }
    for (std::size_t b = 0; b < num_bins; ++b)
      offsets_[b + 1] += offsets_[b];
    cell_ids_.resize(num_cells);
    auto next = offsets_;
    for (std::size_t id = 0; id < num_cells; ++id)
      cell_ids_[next[cell_bin[id]]++] = static_cast<std::uint32_t>(id);
  }

  /// Calls f(local_id) for each local cell whose centroid lies in a bin that overlaps [lo, hi].
  template <typename F>
  void ForEachCandidate(const Vector3& lo, const Vector3& hi, const F& f) const
  {
    if (cell_ids_.empty())
      return;
    for (auto i = BinIndex(0, lo.x); i <= BinIndex(0, hi.x); ++i)
      for (auto j = BinIndex(1, lo.y); j <= BinIndex(1, hi.y); ++j)
        for (auto k = BinIndex(2, lo.z); k <= BinIndex(2, hi.z); ++k)
        {
          const auto b = (i * num_bins_[1] + j) * num_bins_[2] + k;
          for (auto n = offsets_[b]; n < offsets_[b + 1]; ++n)
            f(cell_ids_[n]);
        }
  }

private:
  /**
   * Bin of coordinate x in direction d, clamped to the grid. It is non-decreasing in x, so a
   * centroid between lo and hi falls in a bin between those of lo and hi. Infinite x is allowed.
   */
  std::size_t BinIndex(std::size_t d, double x) const
  {
    if (num_bins_[d] == 1)
      return 0;
    const double t = (x - origin_[d]) / bin_width_[d];
    if (not(t > 0.0))
      return 0;
    if (t >= static_cast<double>(num_bins_[d] - 1))
      return num_bins_[d] - 1;
    return static_cast<std::size_t>(t);
  }

  std::array<double, 3> origin_{};
  std::array<double, 3> bin_width_{};
  std::array<std::size_t, 3> num_bins_{1, 1, 1};
  std::vector<std::size_t> offsets_;
  std::vector<std::uint32_t> cell_ids_;
};

/**
 * Local IDs, in increasing order, of the local cells whose centroid is inside log_vol and whose
 * block ID is in block_ids (any block when block_ids is empty).
 */
std::vector<std::uint32_t>
GetLogicalVolumeCellIDs(const Mesh& grid,
                        const CentroidBins& bins,
                        const LogicalVolume& log_vol,
                        const std::vector<int>& block_ids)
{
  std::vector<std::uint32_t> cell_ids;
  const auto add_if_selected = [&](const std::uint32_t id)
  {
    const auto& cell = grid.GetLocalCell(id);
    if (log_vol.Inside(cell.centroid) and
        (block_ids.empty() or
         std::find(block_ids.begin(), block_ids.end(), cell.block_id) != block_ids.end()))
      cell_ids.push_back(id);
  };

  Vector3 lo, hi;
  if (log_vol.GetBoundingBox(lo, hi))
  {
    bins.ForEachCandidate(lo, hi, add_if_selected);
    std::sort(cell_ids.begin(), cell_ids.end());
  }
  else
  {
    for (const auto& cell : grid.GetLocalCells())
      add_if_selected(cell->local_id);
  }
  return cell_ids;
}

} // namespace

InputParameters
VolumePostprocessor::GetInputParameters()
{
  InputParameters params;
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "A handle to an existing LBS problem.");
  params.AddOptionalParameterArray(
    "block_ids",
    std::vector<int>{},
    "Block restriction for the postprocessor. Empty/unspecified means no block restriction.");
  params.AddOptionalParameterArray("logical_volumes",
                                   std::vector<std::shared_ptr<LogicalVolume>>{},
                                   "Logical volume to restrict the computation to.");
  params.AddOptionalParameter<std::string>(
    "value_type", "integral", "Type of value to compute: 'integral', 'max', 'min', or 'avg'");
  params.AddOptionalParameter(
    "group", 0, "Single group to compute (mutually exclusive with groupset).");
  params.AddOptionalParameter(
    "groupset", 0, "Single groupset to compute (mutually exclusive with group).");

  params.AddOptionalParameter("multiplier", 1., "Constant multiplier.");
  params.AddOptionalParameterArray(
    "group_multipliers", std::vector<double>{}, "Group multipliers, one per group.");
  params.AddOptionalParameter("xs_multiplier", "", "Cross-sections multiplier.");

  params.ConstrainParameterRange("value_type",
                                 AllowableRangeList::New({"integral", "max", "min", "avg"}));

  return params;
}

std::shared_ptr<VolumePostprocessor>
VolumePostprocessor::Create(const ParameterBlock& params)
{
  return CreateObject<VolumePostprocessor>("lbs::VolumePostprocessor", params);
}

VolumePostprocessor::VolumePostprocessor(const InputParameters& params)
  : lbs_problem_(params.GetSharedPtrParam<Problem, LBSProblem>("problem")),
    block_ids_(params.GetParamVectorValue<int>("block_ids")),
    logical_volumes_(params.GetParamVectorValue<std::shared_ptr<LogicalVolume>>("logical_volumes")),
    selected_group_(params.IsParameterValid("group")
                      ? std::make_optional(params.GetParamValue<unsigned int>("group"))
                      : std::nullopt),
    selected_groupset_(params.IsParameterValid("groupset")
                         ? std::make_optional(params.GetParamValue<unsigned int>("groupset"))
                         : std::nullopt),
    const_multiplier_(params.IsParameterValid("multiplier")
                        ? std::make_optional(params.GetParamValue<double>("multiplier"))
                        : std::nullopt),
    group_multipliers_(
      params.IsParameterValid("group_multipliers")
        ? std::make_optional(params.GetParamVectorValue<double>("group_multipliers"))
        : std::nullopt),
    xs_multiplier_(params.IsParameterValid("xs_multiplier")
                     ? std::make_optional(params.GetParamValue<std::string>("xs_multiplier"))
                     : std::nullopt)
{
  if (selected_group_.has_value() && selected_groupset_.has_value())
    throw std::invalid_argument("'group' and 'groupset' cannot be specified together");

  if (selected_group_.has_value() && selected_group_.value() >= lbs_problem_->GetNumGroups())
    throw std::invalid_argument("'group' must be less than " +
                                std::to_string(lbs_problem_->GetNumGroups()));

  if (selected_groupset_.has_value() &&
      selected_groupset_.value() >= lbs_problem_->GetGroupsets().size())
    throw std::invalid_argument("'groupset' must be less than " +
                                std::to_string(lbs_problem_->GetGroupsets().size()));

  const auto value_type_str = params.GetParamValue<std::string>("value_type");
  if (value_type_str == "max")
    value_type_ = ValueType::MAX;
  else if (value_type_str == "min")
    value_type_ = ValueType::MIN;
  else if (value_type_str == "integral")
    value_type_ = ValueType::INTEGRAL;
  else if (value_type_str == "avg")
    value_type_ = ValueType::AVERAGE;
  else
    throw std::invalid_argument("'value_type' can be only 'min', 'max', 'integral', or 'avg'");

  CreateSpatialRestriction();
  CreateEnergyRestriction();
  CreateMultipliers();

  auto n_lvs = std::max(static_cast<std::size_t>(1), logical_volumes_.size());
  values_.resize({n_lvs, groups_.size()});
}

void
VolumePostprocessor::CreateSpatialRestriction()
{
  const auto& grid = lbs_problem_->GetGrid();

  if (logical_volumes_.empty())
  {
    cell_local_ids_.resize(1);
    std::vector<std::uint32_t> cell_ids;
    if (block_ids_.empty())
    {
      for (const auto& cell : grid->GetLocalCells())
        cell_ids.push_back(cell->local_id);
    }
    else
    {
      for (const auto& cell : grid->GetLocalCells())
      {
        if (std::find(block_ids_.begin(), block_ids_.end(), cell->block_id) != block_ids_.end())
          cell_ids.push_back(cell->local_id);
      }
    }
    cell_local_ids_[0] = cell_ids;
  }
  else
  {
    const CentroidBins bins(*grid);
    cell_local_ids_.resize(logical_volumes_.size());
    for (std::size_t i = 0; i < logical_volumes_.size(); ++i)
      cell_local_ids_[i] = GetLogicalVolumeCellIDs(*grid, bins, *logical_volumes_[i], block_ids_);
  }

  // An empty selection would silently report zero (or +/- infinity for min/max).
  std::vector<std::size_t> local_counts(cell_local_ids_.size());
  for (std::size_t i = 0; i < cell_local_ids_.size(); ++i)
    local_counts[i] = cell_local_ids_[i].size();
  std::vector<std::size_t> counts(local_counts.size());
  mpi_comm.all_reduce(local_counts, counts, mpi::op::sum<std::size_t>());
  for (std::size_t i = 0; i < counts.size(); ++i)
    OpenSnInvalidArgumentIf(counts[i] == 0,
                            "VolumePostprocessor: " +
                              (logical_volumes_.empty() ? std::string("the block restriction")
                                                        : "logical volume " + std::to_string(i)) +
                              " selects no mesh cells.");
}

void
VolumePostprocessor::CreateMultipliers()
{
  auto n_groups = lbs_problem_->GetNumGroups();

  if (not const_multiplier_.has_value() and not group_multipliers_.has_value() and
      not xs_multiplier_.has_value())
  {
    multipliers_ = std::vector<double>(n_groups, 1.);
  }
  else if (const_multiplier_.has_value() and not group_multipliers_.has_value() and
           not xs_multiplier_.has_value())
  {
    multipliers_ = std::vector<double>(n_groups, const_multiplier_.value());
  }
  else if (not const_multiplier_.has_value() and group_multipliers_.has_value() and
           not xs_multiplier_.has_value())
  {
    OpenSnInvalidArgumentIf(group_multipliers_.value().size() != n_groups,
                            "Must provide one multiplier per group");
    multipliers_ = group_multipliers_.value();
  }
  else if (not const_multiplier_.has_value() and not group_multipliers_.has_value() and
           xs_multiplier_.has_value())
  {
    // do nothing, the actual values will be pulled from xs object during the post-processor
    // computation
  }
  else
    throw std::invalid_argument(
      "Specify at most one of 'multiplier', 'group_multipliers', and 'xs_multiplier'.");
}

void
VolumePostprocessor::CreateEnergyRestriction()
{
  if (selected_group_.has_value())
  {
    groups_.push_back(selected_group_.value());
  }
  else if (selected_groupset_.has_value())
  {
    const auto& groupset = lbs_problem_->GetGroupsets()[selected_groupset_.value()];
    for (unsigned int g = groupset.first_group; g <= groupset.last_group; ++g)
      groups_.push_back(g);
  }
  else
  {
    for (unsigned int g = 0; g < lbs_problem_->GetNumGroups(); ++g)
      groups_.push_back(g);
  }
}

void
VolumePostprocessor::Execute()
{
  if (xs_multiplier_.has_value())
  {
    // The cross sections can change after construction, so check them on each execution. Every
    // rank holds the same block-to-XS map, so every rank reaches the same result.
    const auto num_groups = lbs_problem_->GetNumGroups();
    bool any_block_has_xs = false;
    for (const auto& [block_id, xs] : lbs_problem_->GetBlockID2XSMap())
    {
      const std::vector<double>* coeffs = nullptr;
      try
      {
        coeffs = xs->GetByName(xs_multiplier_.value());
      }
      catch (const std::logic_error& error)
      {
        // An unknown name is invalid input.
        throw std::invalid_argument("VolumePostprocessor: " + std::string(error.what()));
      }
      OpenSnInvalidArgumentIf(coeffs != nullptr and coeffs->size() != num_groups,
                              "VolumePostprocessor: cross section \"" + xs_multiplier_.value() +
                                "\" on block " + std::to_string(block_id) +
                                " does not have one value per group.");
      any_block_has_xs = any_block_has_xs or coeffs != nullptr;
    }
    // A material without the cross section contributes zero, but a problem without it would
    // silently report zero.
    OpenSnInvalidArgumentIf(not any_block_has_xs,
                            "VolumePostprocessor: no material in the problem has cross section \"" +
                              xs_multiplier_.value() + "\".");
    zero_coefficients_.assign(num_groups, 0.0);
  }

  // Compute the local values of every volume, then reduce them all at once.
  const auto num_volumes = cell_local_ids_.size();
  const auto num_groups = groups_.size();
  const auto num_values = num_volumes * num_groups;
  std::vector<double> local;
  std::vector<double> global;
  switch (value_type_)
  {
    case ValueType::INTEGRAL:
    case ValueType::AVERAGE:
    {
      // For an average, the volumes follow the integrals.
      const bool average = value_type_ == ValueType::AVERAGE;
      local.assign(num_values + (average ? num_volumes : 0), 0.0);
      for (std::size_t i = 0; i < num_volumes; ++i)
        AccumulateLocalIntegral(
          cell_local_ids_[i], &local[i * num_groups], average ? &local[num_values + i] : nullptr);
      global.resize(local.size());
      mpi_comm.all_reduce(
        local.data(), static_cast<int>(local.size()), global.data(), mpi::op::sum<double>());
      for (std::size_t i = 0; i < num_volumes; ++i)
      {
        const double volume = average ? global[num_values + i] : 1.0;
        for (std::size_t k = 0; k < num_groups; ++k)
          values_(i, k) = volume > 0.0 ? global[i * num_groups + k] / volume : 0.0;
      }
      break;
    }
    case ValueType::MAX:
    case ValueType::MIN:
    {
      const bool is_max = value_type_ == ValueType::MAX;
      local.resize(num_values);
      for (std::size_t i = 0; i < num_volumes; ++i)
        ComputeLocalExtremum(cell_local_ids_[i], is_max, &local[i * num_groups]);
      global.resize(num_values);
      if (is_max)
        mpi_comm.all_reduce(
          local.data(), static_cast<int>(num_values), global.data(), mpi::op::max<double>());
      else
        mpi_comm.all_reduce(
          local.data(), static_cast<int>(num_values), global.data(), mpi::op::min<double>());
      for (std::size_t i = 0; i < num_volumes; ++i)
        for (std::size_t k = 0; k < num_groups; ++k)
          values_(i, k) = global[i * num_groups + k];
      break;
    }
  }
}

const std::vector<double>&
VolumePostprocessor::GetCoefficients(const Cell& cell)
{
  if (xs_multiplier_.has_value())
  {
    const auto blk_id = cell.block_id;
    const auto& xs = lbs_problem_->GetBlockID2XSMap().at(blk_id);
    // A material without this cross section (e.g. sigma_f of a non-fissile one) contributes zero.
    const auto* coeffs = xs->GetByName(xs_multiplier_.value());
    return coeffs != nullptr ? *coeffs : zero_coefficients_;
  }
  else
  {
    return multipliers_;
  }
}

void
VolumePostprocessor::AccumulateLocalIntegral(const std::vector<uint32_t>& cell_local_ids,
                                             double* integral,
                                             double* volume)
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto& phi = lbs_problem_->GetPhiNewLocal();
  auto coord = sdm.GetSpatialWeightingFunction();

  for (const auto cell_id : cell_local_ids)
  {
    const auto& cell = grid->GetLocalCell(cell_id);
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto num_nodes = cell_mapping.GetNumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();
    const auto& coeffs = GetCoefficients(cell);

    if (volume != nullptr)
      for (const std::size_t qp : fe_vol_data.GetQuadraturePointIndices())
        *volume += coord(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.JxW(qp);

    for (std::size_t k = 0; k < groups_.size(); ++k)
    {
      std::vector<double> nodal_value(num_nodes, 0.0);
      for (std::size_t i = 0; i < num_nodes; ++i)
      {
        const auto imap = sdm.MapDOFLocal(cell, i, uk_man, 0, groups_[k]);
        nodal_value[i] = phi[imap];
      }

      for (const std::size_t qp : fe_vol_data.GetQuadraturePointIndices())
      {
        double phi_h = 0.0;
        for (std::size_t j = 0; j < num_nodes; ++j)
          phi_h += fe_vol_data.ShapeValue(j, qp) * nodal_value[j];

        integral[k] +=
          coeffs[groups_[k]] * phi_h * coord(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.JxW(qp);
      }
    }
  }
}

void
VolumePostprocessor::ComputeLocalExtremum(const std::vector<uint32_t>& cell_local_ids,
                                          const bool is_max,
                                          double* extremum)
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto& phi = lbs_problem_->GetPhiNewLocal();

  const double init =
    is_max ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
  std::fill(extremum, extremum + groups_.size(), init);
  for (const auto cell_id : cell_local_ids)
  {
    const auto& cell = grid->GetLocalCell(cell_id);
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto num_nodes = cell_mapping.GetNumNodes();
    const auto& coeffs = GetCoefficients(cell);

    for (std::size_t k = 0; k < groups_.size(); ++k)
    {
      for (std::size_t i = 0; i < num_nodes; ++i)
      {
        const auto imap = sdm.MapDOFLocal(cell, i, uk_man, 0, groups_[k]);
        const double value = coeffs[groups_[k]] * phi[imap];
        extremum[k] = is_max ? std::max(extremum[k], value) : std::min(extremum[k], value);
      }
    }
  }
}

const NDArray<double, 2>&
VolumePostprocessor::GetValue() const
{
  return values_;
}

} // namespace opensn
