// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/error.h"

namespace opensn
{
void
SweepChunk::Sweep(AngleSet& angle_set)
{
  OpenSnLogicalError("SweepChunk::Sweep must be overridden by derived sweep chunk.");
}

void
SweepChunk::SetAngleSet(AngleSet& angle_set)
{
  OpenSnLogicalError("SweepChunk::SetAngleSet is not implemented for this sweep chunk.");
}

void
SweepChunk::SetCell(Cell const* cell_ptr, AngleSet& angle_set)
{
  OpenSnLogicalError("SweepChunk::SetCell is not implemented for this sweep chunk.");
}

void
SweepChunk::ZeroDestinationPhi()
{
  const auto gsi = groupset_.first_group;
  const auto gss = groupset_.GetNumGroups();

  for (const auto& cell : grid_->local_cells)
  {
    const auto& transport_view = cell_transport_views_[cell.local_id];
    const auto num_nodes = static_cast<size_t>(transport_view.GetNumNodes());
    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (unsigned int m = 0; m < num_moments_; ++m)
      {
        const auto mapping = transport_view.MapDOF(i, m, gsi);
        for (unsigned int g = 0; g < gss; ++g)
        {
          destination_phi_[mapping + g] = 0.0;
        } // for g
      } // for moment
    } // for dof
  } // for cell
}

} // namespace opensn
