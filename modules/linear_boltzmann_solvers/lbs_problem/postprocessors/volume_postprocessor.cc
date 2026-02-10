// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/volume_postprocessor.h"
#include "framework/object_factory.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/runtime.h"
#include <limits>
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
  params.AddOptionalParameter<std::shared_ptr<LogicalVolume>>(
    "logical_volume", nullptr, "Logical volume to restrict the computation to.");
  params.AddOptionalParameter<std::string>(
    "value_type", "integral", "Type of value to compute: 'integral', 'max', 'min', or 'avg'");
  params.AddOptionalParameter(
    "group", -1, "Single group to compute (mutually exclusive with groupset).");
  params.AddOptionalParameter(
    "groupset", -1, "Single groupset to compute (mutually exclusive with group).");
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
    logical_volume_(params.GetParamValue<std::shared_ptr<LogicalVolume>>("logical_volume")),
    selected_group_(params.IsParameterValid("group")
                      ? std::make_optional(params.GetParamValue<unsigned int>("group"))
                      : std::nullopt),
    selected_groupset_(params.IsParameterValid("groupset")
                         ? std::make_optional(params.GetParamValue<unsigned int>("groupset"))
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
  values_.resize(groups_.size());
}

void
VolumePostprocessor::CreateSpatialRestriction()
{
  const auto& grid = lbs_problem_->GetGrid();

  // filter on logical volumes
  std::vector<std::uint32_t> cell_ids;
  if (logical_volume_ != nullptr)
  {
    for (const auto& cell : grid->local_cells)
      if (logical_volume_->Inside(cell.centroid))
        cell_ids.push_back(cell.local_id);
  }
  else
  {
    for (const auto& cell : grid->local_cells)
      cell_ids.push_back(cell.local_id);
  }

  // apply block restriction
  if (block_ids_.empty())
  {
    cell_local_ids_.assign(cell_ids.begin(), cell_ids.end());
  }
  else
  {
    for (const auto& id : cell_ids)
    {
      auto block_id = grid->local_cells[id].block_id;
      if (std::find(block_ids_.begin(), block_ids_.end(), block_id) != block_ids_.end())
        cell_local_ids_.push_back(id);
    }
  }

  // TODO: if `cell_local_ids_` is empty, warn or stop
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
    for (unsigned int g = groupset.first_group; g < groupset.last_group; ++g)
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
      ComputeIntegral();
      break;
    case ValueType::MAX:
      ComputeMax();
      break;
    case ValueType::MIN:
      ComputeMin();
      break;
    case ValueType::AVERAGE:
      ComputeVolumeWeightedAverage();
      break;
  }
}

void
VolumePostprocessor::ComputeIntegral()
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto phi = lbs_problem_->GetPhiNewLocal();
  auto coord = sdm.GetSpatialWeightingFunction();

  std::vector<double> local_integral(groups_.size(), 0.0);
  for (const auto cell_local_id : cell_local_ids_)
  {
    const auto& cell = grid->local_cells[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto num_nodes = cell_mapping.GetNumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

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

        local_integral[k] += phi_h * coord(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.JxW(qp);
      }
    }
  }

  std::vector<double> global_integral(groups_.size(), 0.0);
  for (std::size_t i = 0; i < local_integral.size(); ++i)
    mpi_comm.all_reduce(local_integral[i], global_integral[i], mpi::op::sum<double>());

  for (std::size_t i = 0; i < global_integral.size(); ++i)
    values_[i] = global_integral[i];
}

void
VolumePostprocessor::ComputeMax()
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto phi = lbs_problem_->GetPhiNewLocal();

  std::vector<double> local_max(groups_.size(), -std::numeric_limits<double>::infinity());
  for (const auto cell_local_id : cell_local_ids_)
  {
    const auto& cell = grid->local_cells[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto num_nodes = cell_mapping.GetNumNodes();

    for (std::size_t k = 0; k < groups_.size(); ++k)
    {
      for (std::size_t i = 0; i < num_nodes; ++i)
      {
        const auto imap = sdm.MapDOFLocal(cell, i, uk_man, 0, groups_[k]);
        local_max[k] = std::max(local_max[k], phi[imap]);
      }
    }
  }

  std::vector<double> global_max(groups_.size(), -std::numeric_limits<double>::infinity());
  for (std::size_t i = 0; i < local_max.size(); ++i)
    mpi_comm.all_reduce(local_max[i], global_max[i], mpi::op::max<double>());

  for (std::size_t i = 0; i < global_max.size(); ++i)
    values_[i] = global_max[i];
}

void
VolumePostprocessor::ComputeMin()
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto phi = lbs_problem_->GetPhiNewLocal();

  std::vector<double> local_min(groups_.size(), std::numeric_limits<double>::infinity());
  for (const auto cell_local_id : cell_local_ids_)
  {
    const auto& cell = grid->local_cells[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto num_nodes = cell_mapping.GetNumNodes();

    for (std::size_t k = 0; k < groups_.size(); ++k)
    {
      for (std::size_t i = 0; i < num_nodes; ++i)
      {
        const auto imap = sdm.MapDOFLocal(cell, i, uk_man, 0, groups_[k]);
        local_min[k] = std::min(local_min[k], phi[imap]);
      }
    }
  }

  std::vector<double> global_min(groups_.size(), std::numeric_limits<double>::infinity());
  for (std::size_t i = 0; i < local_min.size(); ++i)
    mpi_comm.all_reduce(local_min[i], global_min[i], mpi::op::min<double>());

  for (std::size_t i = 0; i < global_min.size(); ++i)
    values_[i] = global_min[i];
}

void
VolumePostprocessor::ComputeVolumeWeightedAverage()
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto phi = lbs_problem_->GetPhiNewLocal();
  auto coord = sdm.GetSpatialWeightingFunction();

  std::vector<double> local_weighted_integral(groups_.size(), 0.0);
  double local_weighted_volume = 0.0;

  for (const auto cell_local_id : cell_local_ids_)
  {
    const auto& cell = grid->local_cells[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto num_nodes = cell_mapping.GetNumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

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

        const auto weight = coord(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.JxW(qp);
        local_weighted_integral[k] += phi_h * weight;
      }
    }

    for (const std::size_t qp : fe_vol_data.GetQuadraturePointIndices())
      local_weighted_volume += coord(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.JxW(qp);
  }

  std::vector<double> global_weighted_integral(groups_.size(), 0.0);
  double global_weighted_volume = 0.0;

  for (std::size_t i = 0; i < local_weighted_integral.size(); ++i)
    mpi_comm.all_reduce(
      local_weighted_integral[i], global_weighted_integral[i], mpi::op::sum<double>());
  mpi_comm.all_reduce(local_weighted_volume, global_weighted_volume, mpi::op::sum<double>());

  for (std::size_t i = 0; i < global_weighted_integral.size(); ++i)
  {
    if (global_weighted_volume > 0.0)
      values_[i] = global_weighted_integral[i] / global_weighted_volume;
    else
      values_[i] = 0.0;
  }
}

std::vector<double>
VolumePostprocessor::GetValue() const
{
  return values_;
}

} // namespace opensn
