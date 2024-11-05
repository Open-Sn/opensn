// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/interpolation/ffinter_point.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{

void
FieldFunctionInterpolationPoint::Initialize()
{
  log.Log0Verbose1() << "Initializing point interpolator.";

  // Check for empty FF-list
  if (field_functions_.empty())
    throw std::logic_error("Unassigned field function in point field function interpolator.");

  const auto& grid = field_functions_.front()->GetSpatialDiscretization().Grid();
  std::vector<uint64_t> cells_potentially_owning_point;
  for (const auto& cell : grid.local_cells)
  {
    const auto& vcc = cell.centroid;
    const auto& poi = point_of_interest_;
    const auto nudged_point = poi + 1.0e-6 * (vcc - poi);
    if (grid.CheckPointInsideCell(cell, nudged_point))
      cells_potentially_owning_point.push_back(cell.global_id);
  }

  std::vector<uint64_t> recvbuf;
  mpi_comm.all_gather(cells_potentially_owning_point, recvbuf);

  if (recvbuf.empty())
  {
    throw std::logic_error("FieldFunctionInterpolationPoint::Initialize: No cell identified "
                           "containing the specified point.");
  }

  uint64_t owning_cell_gid = recvbuf.front();
  for (const uint64_t gid : recvbuf)
    owning_cell_gid = std::min(owning_cell_gid, gid);
  locally_owned_ = false;
  for (const uint64_t gid : cells_potentially_owning_point)
  {
    if (gid == owning_cell_gid)
    {
      locally_owned_ = true;
      owning_cell_gid_ = owning_cell_gid;
      break;
    }
  }
}

void
FieldFunctionInterpolationPoint::Execute()
{
  if (not locally_owned_)
    return;

  const auto& ref_ff = *field_functions_.front();
  const auto& sdm = ref_ff.GetSpatialDiscretization();
  const auto& grid = sdm.Grid();

  const auto& uk_man = ref_ff.GetUnknownManager();
  const auto uid = 0;
  const auto cid = ref_component_;

  const auto field_data = ref_ff.GhostedFieldVector();

  const auto& cell = grid.cells[owning_cell_gid_];
  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();

  std::vector<double> node_dof_values(num_nodes, 0.0);
  for (size_t i = 0; i < num_nodes; ++i)
  {
    const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);
    node_dof_values[i] = field_data[imap];
  }

  Vector<double> shape_values(num_nodes, 0.0);
  cell_mapping.ShapeValues(point_of_interest_, shape_values);
  point_value_ = 0.0;
  for (size_t i = 0; i < num_nodes; ++i)
    point_value_ += node_dof_values[i] * shape_values(i);
}

double
FieldFunctionInterpolationPoint::PointValue() const
{
  double global_point_value;
  mpi_comm.all_reduce(point_value_, global_point_value, mpi::op::sum<double>());
  return global_point_value;
}

} // namespace opensn
