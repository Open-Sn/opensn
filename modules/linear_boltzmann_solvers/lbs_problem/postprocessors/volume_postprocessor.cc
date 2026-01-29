// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/volume_postprocessor.h"
#include "framework/object_factory.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/runtime.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, VolumePostprocessor);

InputParameters
VolumePostprocessor::GetInputParameters()
{
  InputParameters params;
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "A handle to an existing LBS problem.");
  params.AddOptionalParameter<std::shared_ptr<LogicalVolume>>("logical_volume",
                                                               nullptr,
                                                               "Logical volume to restrict the computation to.");
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
    logical_volume_(params.GetParamValue<std::shared_ptr<LogicalVolume>>("logical_volume"))
{
}

void
VolumePostprocessor::Initialize()
{
  const auto& grid = lbs_problem_->GetGrid();

  if (logical_volume_ != nullptr)
  {
    for (const auto& cell : grid->local_cells)
      if (logical_volume_->Inside(cell.centroid))
        cell_local_ids_.push_back(cell.local_id);
  }
  else
  {
    for (const auto& cell : grid->local_cells)
      cell_local_ids_.push_back(cell.local_id);
  }

  // TODO: do a single groupset, or a single group
  for (unsigned int g = 0; g < lbs_problem_->GetNumGroups(); ++g)
    groups_.push_back(g);

  values_.resize(groups_.size());
}

void
VolumePostprocessor::Execute()
{
  const auto& sdm = lbs_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();

  const auto& uk_man = lbs_problem_->GetUnknownManager();
  const auto uid = 0;
  const auto cid = 0;

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
        // phi_h = sum_j b_j phi_j
        double phi_h = 0.0;
        for (std::size_t j = 0; j < num_nodes; ++j)
          phi_h += fe_vol_data.ShapeValue(j, qp) * nodal_value[j];

        local_integral[k] += phi_h * coord(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.JxW(qp);
      }
    }
  }

  std::vector<double> global_integral(groups_.size(), 0.0);
  for (std::size_t i = 0; i < local_integral.size(); i++)
    mpi_comm.all_reduce(local_integral[i], global_integral[i], mpi::op::sum<double>());

  for (std::size_t i = 0; i < global_integral.size(); i++)
  {
    std::cerr << "pps = " << global_integral[i] << std::endl;
    values_[i] = global_integral[i];
  }
}

std::vector<double>
VolumePostprocessor::GetValue() const
{
  return values_;
}

} // namespace opensn
