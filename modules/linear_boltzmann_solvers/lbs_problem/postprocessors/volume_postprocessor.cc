// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/volume_postprocessor.h"
#include "framework/object_factory.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep_parallel_for.h"
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, VolumePostprocessor);

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
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<VolumePostprocessor>("lbs::VolumePostprocessor", params);
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
      for (const auto& cell : grid->local_cells)
        cell_ids.push_back(cell.local_id);
    }
    else
    {
      for (const auto& cell : grid->local_cells)
      {
        if (std::find(block_ids_.begin(), block_ids_.end(), cell.block_id) != block_ids_.end())
          cell_ids.push_back(cell.local_id);
      }
    }
    cell_local_ids_[0] = cell_ids;
  }
  else
  {
    cell_local_ids_.resize(logical_volumes_.size());
    for (unsigned int i = 0; i < logical_volumes_.size(); ++i)
      cell_local_ids_[i] = GetLogicalVolumeCellIDs(logical_volumes_[i]);
  }

  OpenSnLogicalErrorIf(cell_local_ids_.empty(),
                       "No mesh cells were selected by a VolumePostprocessor.");
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
    throw std::logic_error("Can specify either 'multiplier' or 'group_multipliers', not both.");
}

std::vector<std::uint32_t>
VolumePostprocessor::GetLogicalVolumeCellIDs(std::shared_ptr<LogicalVolume> log_vol)
{
  const auto& grid = lbs_problem_->GetGrid();

  // filter on logical volumes
  std::vector<std::uint32_t> cell_ids;
  for (const auto& cell : grid->local_cells)
    if (log_vol->Inside(cell.centroid))
      cell_ids.push_back(cell.local_id);

  std::vector<std::uint32_t> final_cell_ids;
  // apply block restriction
  if (block_ids_.empty())
  {
    final_cell_ids.assign(cell_ids.begin(), cell_ids.end());
  }
  else
  {
    for (const auto& id : cell_ids)
    {
      auto block_id = grid->local_cells[id].block_id;
      if (std::find(block_ids_.begin(), block_ids_.end(), block_id) != block_ids_.end())
        final_cell_ids.push_back(id);
    }
  }

  return final_cell_ids;
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
  switch (value_type_)
  {
    case ValueType::INTEGRAL:
      ExecuteReduction<&VolumePostprocessor::ComputeIntegral>(
        0.0, [](double a, double b) { return a + b; }, mpi::op::sum<double>());
      break;
    case ValueType::MAX:
      ExecuteReduction<&VolumePostprocessor::ComputeMax>(
        -std::numeric_limits<double>::infinity(),
        [](double a, double b) { return std::max(a, b); },
        mpi::op::max<double>());
      break;
    case ValueType::MIN:
      ExecuteReduction<&VolumePostprocessor::ComputeMin>(
        std::numeric_limits<double>::infinity(),
        [](double a, double b) { return std::min(a, b); },
        mpi::op::min<double>());
      break;
    case ValueType::AVERAGE:
      ExecuteWeightedAverage();
      break;
  }
}

template <auto ComputeLocal, typename Combine, typename Op>
void
VolumePostprocessor::ExecuteReduction(double identity, Combine combine, Op op)
{
  const auto n_lvs = cell_local_ids_.size();
  const auto n_groups = groups_.size();

  std::vector<std::pair<std::size_t, std::uint32_t>> work;
  for (std::size_t i = 0; i < n_lvs; ++i)
    for (const auto cell_id : cell_local_ids_[i])
      work.emplace_back(i, cell_id);

  const auto n_threads =
    std::min<std::size_t>(std::max(1U, opensn_num_threads), std::max<std::size_t>(work.size(), 1));

  std::vector<std::vector<double>> thread_local_all(
    n_threads, std::vector<double>(n_lvs * n_groups, identity));

  const auto& grid = lbs_problem_->GetSpatialDiscretization().GetGrid();
  ParallelFor(work.size(),
              n_threads,
              [&](std::size_t w)
              {
                const auto [lv_index, cell_id] = work[w];
                auto& buffer = thread_local_all[w % n_threads];
                std::span<double> row(buffer.data() + lv_index * n_groups, n_groups);
                (this->*ComputeLocal)(row, grid->local_cells[cell_id]);
              });

  std::vector<double> local_all(n_lvs * n_groups, identity);
  for (const auto& buffer : thread_local_all)
    for (std::size_t idx = 0; idx < local_all.size(); ++idx)
      local_all[idx] = combine(local_all[idx], buffer[idx]);

  std::vector<double> global_all(n_lvs * n_groups, identity);
  if (not local_all.empty())
    mpi_comm.all_reduce(local_all, global_all, op);

  for (std::size_t i = 0; i < n_lvs; ++i)
    MemCopyRow(
      values_,
      i,
      std::vector<double>(global_all.begin() +
                            static_cast<std::vector<double>::difference_type>(i * n_groups),
                          global_all.begin() +
                            static_cast<std::vector<double>::difference_type>((i + 1) * n_groups)));
}

void
VolumePostprocessor::ExecuteWeightedAverage()
{
  const auto n_lvs = cell_local_ids_.size();
  const auto n_groups = groups_.size();
  const auto stride = n_groups + 1;

  std::vector<std::pair<std::size_t, std::uint32_t>> work;
  for (std::size_t i = 0; i < n_lvs; ++i)
    for (const auto cell_id : cell_local_ids_[i])
      work.emplace_back(i, cell_id);

  const auto n_threads =
    std::min<std::size_t>(std::max(1U, opensn_num_threads), std::max<std::size_t>(work.size(), 1));

  std::vector<std::vector<double>> thread_local_all(n_threads,
                                                    std::vector<double>(n_lvs * stride, 0.0));

  const auto& grid = lbs_problem_->GetSpatialDiscretization().GetGrid();
  ParallelFor(work.size(),
              n_threads,
              [&](std::size_t w)
              {
                const auto [lv_index, cell_id] = work[w];
                const auto& cell = grid->local_cells[cell_id];
                auto& buffer = thread_local_all[w % n_threads];
                buffer[lv_index * stride] += ComputeWeightedVolume(cell);
                std::span<double> row(buffer.data() + lv_index * stride + 1, n_groups);
                ComputeIntegral(row, cell);
              });

  std::vector<double> local_all(n_lvs * stride, 0.0);
  for (const auto& buffer : thread_local_all)
    for (std::size_t idx = 0; idx < local_all.size(); ++idx)
      local_all[idx] += buffer[idx];

  std::vector<double> global_all(n_lvs * stride, 0.0);
  if (not local_all.empty())
    mpi_comm.all_reduce(local_all, global_all, mpi::op::sum<double>());

  for (std::size_t i = 0; i < n_lvs; ++i)
  {
    const double weighted_volume = global_all[i * stride];
    std::vector<double> values(n_groups, 0.0);
    for (std::size_t k = 0; k < n_groups; ++k)
    {
      const double weighted_integral = global_all[i * stride + 1 + k];
      values[k] = weighted_volume > 0.0 ? weighted_integral / weighted_volume : 0.0;
    }
    MemCopyRow(values_, i, values);
  }
}

const std::vector<double>&
VolumePostprocessor::GetCoefficients(const Cell& cell)
{
  if (xs_multiplier_.has_value())
  {
    const auto blk_id = cell.block_id;
    const auto xs_map = lbs_problem_->GetBlockID2XSMap();
    const auto& xs = xs_map.at(blk_id);
    const auto* coeffs = xs->GetByName(xs_multiplier_.value());
    if (coeffs == nullptr)
      throw std::runtime_error("Requested cross-section name does not exist.");
    return *coeffs;
  }
  else
  {
    return multipliers_;
  }
}

void
VolumePostprocessor::ComputeIntegral(std::span<double> row, const Cell& cell)
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto phi = lbs_problem_->GetPhiNewLocal();
  auto coord = sdm.GetSpatialWeightingFunction();

  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const auto num_nodes = cell_mapping.GetNumNodes();
  const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();
  const auto& coeffs = GetCoefficients(cell);

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

      row[k] += coeffs[groups_[k]] * phi_h * coord(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.JxW(qp);
    }
  }
}

void
VolumePostprocessor::ComputeMax(std::span<double> row, const Cell& cell)
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto phi = lbs_problem_->GetPhiNewLocal();

  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const auto num_nodes = cell_mapping.GetNumNodes();
  const auto& coeffs = GetCoefficients(cell);

  for (std::size_t k = 0; k < groups_.size(); ++k)
    for (std::size_t i = 0; i < num_nodes; ++i)
    {
      const auto imap = sdm.MapDOFLocal(cell, i, uk_man, 0, groups_[k]);
      row[k] = std::max(row[k], coeffs[groups_[k]] * phi[imap]);
    }
}

void
VolumePostprocessor::ComputeMin(std::span<double> row, const Cell& cell)
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto phi = lbs_problem_->GetPhiNewLocal();

  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const auto num_nodes = cell_mapping.GetNumNodes();
  const auto& coeffs = GetCoefficients(cell);

  for (std::size_t k = 0; k < groups_.size(); ++k)
    for (std::size_t i = 0; i < num_nodes; ++i)
    {
      const auto imap = sdm.MapDOFLocal(cell, i, uk_man, 0, groups_[k]);
      row[k] = std::min(row[k], coeffs[groups_[k]] * phi[imap]);
    }
}

double
VolumePostprocessor::ComputeWeightedVolume(const Cell& cell)
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  auto coord = sdm.GetSpatialWeightingFunction();
  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

  double weighted_volume = 0.0;
  for (const std::size_t qp : fe_vol_data.GetQuadraturePointIndices())
    weighted_volume += coord(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.JxW(qp);

  return weighted_volume;
}

const NDArray<double, 2>&
VolumePostprocessor::GetValue() const
{
  return values_;
}

} // namespace opensn
