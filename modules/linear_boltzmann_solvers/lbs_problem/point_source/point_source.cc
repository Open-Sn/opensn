// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/point_source/point_source.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/functions/function.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <numeric>
#include <limits>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, PointSource);

InputParameters
PointSource::GetInputParameters()
{
  InputParameters params;

  params.SetGeneralDescription("A multi-group isotropic point source.");
  params.SetClassName("Point Source");

  params.AddRequiredParameterArray("location", "The (x, y, z) coordinate of the point source.");
  params.AddOptionalParameterArray(
    "strength", std::vector<double>(), "The group-wise point source strength");
  params.AddOptionalParameter<std::shared_ptr<GroupTimeFunction>>(
    "strength_function",
    std::shared_ptr<GroupTimeFunction>{},
    "Function defining group-wise strengths as a function of (group, time).");
  params.AddOptionalParameter("start_time",
                              -std::numeric_limits<double>::infinity(),
                              "Time at which the source becomes active.");
  params.AddOptionalParameter("end_time",
                              std::numeric_limits<double>::infinity(),
                              "Time at which the source becomes inactive.");

  return params;
}

std::shared_ptr<PointSource>
PointSource::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<PointSource>("lbs::PointSource", params);
}

PointSource::PointSource(const InputParameters& params)
  : location_(params.GetParamVectorValue<double>("location")),
    strength_(params.GetParamVectorValue<double>("strength")),
    strength_function_(params.GetSharedPtrParam<GroupTimeFunction>("strength_function", false)),
    start_time_(params.GetParamValue<double>("start_time")),
    end_time_(params.GetParamValue<double>("end_time"))
{
  const bool has_strength = not strength_.empty();
  const bool has_strength_func = static_cast<bool>(strength_function_);
  if ((has_strength ? 1 : 0) + (has_strength_func ? 1 : 0) == 0)
    throw std::invalid_argument("Either a strength vector or strength_function must be provided.");
  if ((has_strength ? 1 : 0) + (has_strength_func ? 1 : 0) > 1)
    throw std::invalid_argument(
      "Specify only one of strength or strength_function for a point source.");
  if (has_strength_func && (start_time_ != -std::numeric_limits<double>::infinity() ||
                            end_time_ != std::numeric_limits<double>::infinity()))
    throw std::invalid_argument("strength_function cannot be used with start_time/end_time. "
                                "Define time dependence in the callback or omit the time bounds.");
  if (not strength_.empty() &&
      std::all_of(strength_.begin(), strength_.end(), [](double x) { return x == 0.0; }))
    log.Log0Warning() << "Point source at " << location_.PrintStr() << " "
                      << "does not have a non-zero source strength.";
}

void
PointSource::Initialize(const LBSProblem& lbs_problem)
{
  if (not strength_function_)
  {
    OpenSnLogicalErrorIf(strength_.size() != lbs_problem.GetNumGroups(),
                         "Incompatible point source strength vector at location " +
                           location_.PrintStr() + ". " + "There are " +
                           std::to_string(lbs_problem.GetNumGroups()) + " energy groups, but " +
                           std::to_string(strength_.size()) + " source strength values.");
  }

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

bool
PointSource::IsActive(double time) const
{
  return time >= start_time_ && time <= end_time_;
}

std::vector<double>
PointSource::GetStrength(const double time, const unsigned int num_groups) const
{
  if (strength_function_)
  {
    std::vector<double> values(num_groups, 0.0);
    for (unsigned int g = 0; g < num_groups; ++g)
      values[g] = (*strength_function_)(g, time);
    return values;
  }
  return strength_;
}

} // namespace opensn
