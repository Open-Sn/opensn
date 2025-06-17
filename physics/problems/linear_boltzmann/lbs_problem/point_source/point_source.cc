// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/problems/linear_boltzmann/lbs_problem/point_source/point_source.h"
#include "physics/problems/linear_boltzmann/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include <numeric>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, PointSource);

InputParameters
PointSource::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("A multi-group isotropic point source.");
  params.SetClassName("Point Source");
  params.SetDocGroup("LBSUtilities");

  params.AddRequiredParameterArray("location", "The (x, y, z) coordinate of the point source.");
  params.AddRequiredParameterArray("strength", "The group-wise point source strength");

  return params;
}

std::shared_ptr<PointSource>
PointSource::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<PointSource>("lbs::PointSource", params);
}

PointSource::PointSource(const InputParameters& params)
  : Object(params),
    location_(params.GetParamVectorValue<double>("location")),
    strength_(params.GetParamVectorValue<double>("strength"))
{
  if (std::all_of(strength_.begin(), strength_.end(), [](double x) { return x == 0.0; }))
    log.Log0Warning() << "Point source at " << location_.PrintStr() << " "
                      << "does not have a non-zero source strength.";
}

void
PointSource::Initialize(const LBSProblem& lbs_problem)
{
  OpenSnLogicalErrorIf(strength_.size() != lbs_problem.GetNumGroups(),
                       "Incompatible point source strength vector at location " +
                         location_.PrintStr() + ". " + "There are " +
                         std::to_string(lbs_problem.GetNumGroups()) + " simulation groups, but " +
                         std::to_string(strength_.size()) + " source strength values.");

  // Get info from solver
  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  const auto& ghost_unit_cell_matrices = lbs_problem.GetUnitGhostCellMatrices();

  // Find local subscribers
  double total_volume = 0.0;
  std::vector<Subscriber> subscribers;
  for (const auto& cell : grid->local_cells)
  {
    if (grid->CheckPointInsideCell(cell, location_))
    {
      const auto& cell_mapping = discretization.GetCellMapping(cell);
      const auto& fe_values = unit_cell_matrices[cell.local_id];

      // Map the point source to the finite element space
      Vector<double> shape_vals;
      cell_mapping.ShapeValues(location_, shape_vals);
      const auto M_inv = Inverse(fe_values.intV_shapeI_shapeJ);
      const auto node_wgts = Mult(M_inv, shape_vals);

      // Increment the total volume
      total_volume += cell.volume;

      // Add to subscribers
      subscribers.push_back(Subscriber{cell.volume, cell.local_id, shape_vals, node_wgts});
    }
  }

  // If the point source lies on a partition boundary, ghost cells must be
  // added to the total volume.
  auto ghost_global_ids = grid->cells.GetGhostGlobalIDs();
  for (uint64_t global_id : ghost_global_ids)
  {
    const auto& nbr_cell = grid->cells[global_id];
    if (grid->CheckPointInsideCell(nbr_cell, location_))
    {
      const auto& fe_values = ghost_unit_cell_matrices.at(nbr_cell.global_id);
      total_volume +=
        std::accumulate(fe_values.intV_shapeI.begin(), fe_values.intV_shapeI.end(), 0.0);
    }
  }

  // Create the actual subscriber list
  subscribers_.clear();
  for (const auto& sub : subscribers)
  {
    subscribers_.push_back(
      {sub.volume_weight / total_volume, sub.cell_local_id, sub.shape_values, sub.node_weights});

    std::stringstream ss;
    ss << "Point source at " << location_.PrintStr() << " assigned to cell "
       << grid->local_cells[sub.cell_local_id].global_id << " with shape values [ ";
    for (const auto& value : sub.shape_values)
      ss << value << " ";
    ss << "] and volume weight " << sub.volume_weight / total_volume;
    log.LogAll() << ss.str();
  }

  size_t num_local_subs = subscribers_.size();
  mpi_comm.all_reduce(num_local_subs, num_global_subscribers_, mpi::op::sum<size_t>());

  log.LogAll() << "Point source has " << num_local_subs << " subscribing cells on processor "
               << mpi_comm.rank();
  log.Log() << "Point source has " << num_global_subscribers_ << " global subscribing cells.";
}

} // namespace opensn
