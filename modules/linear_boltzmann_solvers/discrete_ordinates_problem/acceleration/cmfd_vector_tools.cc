// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_vector_tools.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"

namespace opensn
{

std::vector<double>
CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                       const LBSGroupset& groupset,
                       const CMFDCoarseMesh& coarse_mesh,
                       const std::vector<double>& phi)
{
  return CMFDRestrictScalarFlux(
    do_problem, groupset.first_group, groupset.GetNumGroups(), coarse_mesh, phi);
}

std::vector<double>
CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                       const unsigned int first_group,
                       const unsigned int num_groups,
                       const CMFDCoarseMesh& coarse_mesh,
                       const std::vector<double>& phi)
{
  const auto& transport_views = do_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();

  OpenSnInvalidArgumentIf(phi.size() != do_problem.GetPhiNewLocal().size(),
                          "Input scalar flux vector size does not match the transport problem.");

  std::vector<double> coarse_phi(coarse_mesh.NumLocalCells() * num_groups, 0.0);

  for (const auto& coarse_cell : coarse_mesh.LocalCells())
  {
    std::vector<double> group_integrals(num_groups, 0.0);
    double volume = 0.0;

    for (const uint64_t fine_cell_id : coarse_cell.fine_cell_ids)
    {
      const auto& fine_cell = do_problem.GetGrid()->cells[fine_cell_id];
      OpenSnInvalidArgumentIf(fine_cell.partition_id != opensn::mpi_comm.rank(),
                              "CMFD restriction currently requires local fine cells.");

      const auto& transport_view = transport_views[fine_cell.local_id];
      const auto& cell_matrices = unit_cell_matrices[fine_cell.local_id];
      volume += fine_cell.volume;

      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        const auto node_volume = cell_matrices.intV_shapeI(i);
        const auto phi_map = transport_view.MapDOF(i, 0, first_group);
        for (unsigned int gsg = 0; gsg < num_groups; ++gsg)
          group_integrals[gsg] += phi[phi_map + gsg] * node_volume;
      }
    }

    OpenSnLogicalErrorIf(volume <= 0.0, "CMFD coarse cell has non-positive volume.");

    const auto coarse_offset = coarse_cell.local_id * num_groups;
    for (unsigned int gsg = 0; gsg < num_groups; ++gsg)
      coarse_phi[coarse_offset + gsg] = group_integrals[gsg] / volume;
  }

  return coarse_phi;
}

void
CMFDProlongateScalarFluxCorrection(const DiscreteOrdinatesProblem& do_problem,
                                   const LBSGroupset& groupset,
                                   const CMFDCoarseMesh& coarse_mesh,
                                   const std::vector<double>& coarse_delta_phi,
                                   std::vector<double>& phi)
{
  CMFDProlongateScalarFluxCorrection(
    do_problem, groupset.first_group, groupset.GetNumGroups(), coarse_mesh, coarse_delta_phi, phi);
}

void
CMFDProlongateScalarFluxCorrection(const DiscreteOrdinatesProblem& do_problem,
                                   const unsigned int first_group,
                                   const unsigned int num_groups,
                                   const CMFDCoarseMesh& coarse_mesh,
                                   const std::vector<double>& coarse_delta_phi,
                                   std::vector<double>& phi)
{
  const auto& transport_views = do_problem.GetCellTransportViews();

  OpenSnInvalidArgumentIf(coarse_delta_phi.size() != coarse_mesh.NumLocalCells() * num_groups,
                          "Coarse scalar flux correction vector size mismatch.");
  OpenSnInvalidArgumentIf(phi.size() != do_problem.GetPhiNewLocal().size(),
                          "Output scalar flux vector size does not match the transport problem.");

  for (const auto& coarse_cell : coarse_mesh.LocalCells())
  {
    const auto coarse_offset = coarse_cell.local_id * num_groups;

    for (const uint64_t fine_cell_id : coarse_cell.fine_cell_ids)
    {
      const auto& fine_cell = do_problem.GetGrid()->cells[fine_cell_id];
      OpenSnInvalidArgumentIf(fine_cell.partition_id != opensn::mpi_comm.rank(),
                              "CMFD prolongation currently requires local fine cells.");

      const auto& transport_view = transport_views[fine_cell.local_id];
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        const auto phi_map = transport_view.MapDOF(i, 0, first_group);
        for (unsigned int gsg = 0; gsg < num_groups; ++gsg)
          phi[phi_map + gsg] += coarse_delta_phi[coarse_offset + gsg];
      }
    }
  }
}

} // namespace opensn
