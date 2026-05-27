// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_vector_tools.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <cmath>

namespace opensn
{
namespace
{

unsigned int
NumCoarseGroups(const unsigned int num_groups, const unsigned int group_aggregation_size)
{
  return (num_groups + group_aggregation_size - 1) / group_aggregation_size;
}

unsigned int
CoarseGroupBegin(const unsigned int coarse_group, const unsigned int group_aggregation_size)
{
  return coarse_group * group_aggregation_size;
}

unsigned int
CoarseGroupEnd(const unsigned int coarse_group,
               const unsigned int num_groups,
               const unsigned int group_aggregation_size)
{
  return std::min(num_groups, (coarse_group + 1) * group_aggregation_size);
}

} // namespace

std::vector<double>
CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                       const unsigned int first_group,
                       const unsigned int num_groups,
                       const CMFDCoarseMesh& coarse_mesh,
                       const std::vector<double>& phi)
{
  return CMFDRestrictScalarFlux(do_problem, first_group, num_groups, 1, coarse_mesh, phi);
}

std::vector<double>
CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                       const unsigned int first_group,
                       const unsigned int num_groups,
                       const unsigned int group_aggregation_size,
                       const CMFDCoarseMesh& coarse_mesh,
                       const std::vector<double>& phi)
{
  const auto& transport_views = do_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();

  OpenSnInvalidArgumentIf(phi.size() != do_problem.GetPhiNewLocal().size(),
                          "Input scalar flux vector size does not match the transport problem.");

  OpenSnInvalidArgumentIf(group_aggregation_size == 0,
                          "CMFD group aggregation size must be greater than zero.");

  const unsigned int num_coarse_groups = NumCoarseGroups(num_groups, group_aggregation_size);
  std::vector<double> coarse_phi(coarse_mesh.NumLocalCells() * num_coarse_groups, 0.0);

  for (const auto& coarse_cell : coarse_mesh.LocalCells())
  {
    std::vector<double> group_integrals(num_coarse_groups, 0.0);
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
        for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
          for (unsigned int gsg = CoarseGroupBegin(cg, group_aggregation_size);
               gsg < CoarseGroupEnd(cg, num_groups, group_aggregation_size);
               ++gsg)
            group_integrals[cg] += phi[phi_map + gsg] * node_volume;
      }
    }

    OpenSnLogicalErrorIf(volume <= 0.0, "CMFD coarse cell has non-positive volume.");

    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups;
    for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
      coarse_phi[coarse_offset + cg] = group_integrals[cg] / volume;
  }

  return coarse_phi;
}

void
CMFDProlongateScalarFluxRatio(const DiscreteOrdinatesProblem& do_problem,
                              const unsigned int first_group,
                              const unsigned int num_groups,
                              const unsigned int group_aggregation_size,
                              const CMFDCoarseMesh& coarse_mesh,
                              const std::vector<double>& coarse_phi_old,
                              const std::vector<double>& coarse_phi_new,
                              std::vector<double>& phi)
{
  const auto& transport_views = do_problem.GetCellTransportViews();

  OpenSnInvalidArgumentIf(group_aggregation_size == 0,
                          "CMFD group aggregation size must be greater than zero.");
  const unsigned int num_coarse_groups = NumCoarseGroups(num_groups, group_aggregation_size);
  OpenSnInvalidArgumentIf(coarse_phi_old.size() != coarse_mesh.NumLocalCells() * num_coarse_groups,
                          "Old coarse scalar flux vector size mismatch.");
  OpenSnInvalidArgumentIf(coarse_phi_new.size() != coarse_mesh.NumLocalCells() * num_coarse_groups,
                          "New coarse scalar flux vector size mismatch.");
  OpenSnInvalidArgumentIf(phi.size() != do_problem.GetPhiNewLocal().size(),
                          "Output scalar flux vector size does not match the transport problem.");

  for (const auto& coarse_cell : coarse_mesh.LocalCells())
  {
    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups;

    for (const uint64_t fine_cell_id : coarse_cell.fine_cell_ids)
    {
      const auto& fine_cell = do_problem.GetGrid()->cells[fine_cell_id];
      OpenSnInvalidArgumentIf(fine_cell.partition_id != opensn::mpi_comm.rank(),
                              "CMFD prolongation currently requires local fine cells.");

      const auto& transport_view = transport_views[fine_cell.local_id];
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        const auto phi_map = transport_view.MapDOF(i, 0, first_group);
        for (unsigned int cg = 0; cg < num_coarse_groups; ++cg)
        {
          const double denominator = coarse_phi_old[coarse_offset + cg];
          const double ratio = std::fabs(denominator) > 1.0e-14
                                 ? coarse_phi_new[coarse_offset + cg] / denominator
                                 : 1.0;
          for (unsigned int gsg = CoarseGroupBegin(cg, group_aggregation_size);
               gsg < CoarseGroupEnd(cg, num_groups, group_aggregation_size);
               ++gsg)
            phi[phi_map + gsg] *= ratio;
        }
      }
    }
  }
}

} // namespace opensn
