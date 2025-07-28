// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{

double
ComputeFissionProduction(LBSProblem& lbs_problem, const std::vector<double>& phi)
{
  CALI_CXX_MARK_SCOPE("ComputeFissionProduction");

  auto& groups = lbs_problem.GetGroups();
  auto& grid = lbs_problem.GetGrid();
  auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  auto& options = lbs_problem.GetOptions();

  const int first_grp = groups.front().id;
  const int last_grp = groups.back().id;

  // Loop over local cells
  double local_production = 0.0;
  for (auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& cell_matrices = unit_cell_matrices[cell.local_id];

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& F = xs.GetProductionMatrix();
    const auto& nu_delayed_sigma_f = xs.GetNuDelayedSigmaF();

    if (not xs.IsFissionable())
      continue;

    // Loop over nodes
    const int num_nodes = transport_view.GetNumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.intV_shapeI(i);

      // Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
      {
        const auto& prod = F[g];
        for (size_t gp = 0; gp <= last_grp; ++gp)
          local_production += prod[gp] * phi[uk_map + gp] * IntV_ShapeI;

        if (options.use_precursors)
          for (unsigned int j = 0; j < xs.GetNumPrecursors(); ++j)
            local_production += nu_delayed_sigma_f[g] * phi[uk_map + g] * IntV_ShapeI;
      }
    } // for node
  } // for cell

  // Allreduce global production
  double global_production = 0.0;
  mpi_comm.all_reduce(local_production, global_production, mpi::op::sum<double>());

  return global_production;
}

double
ComputeFissionRate(LBSProblem& lbs_problem, const std::vector<double>& phi)
{
  CALI_CXX_MARK_SCOPE("ComputeFissionRate");

  auto& groups = lbs_problem.GetGroups();
  auto& grid = lbs_problem.GetGrid();
  auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();

  const int first_grp = groups.front().id;
  const int last_grp = groups.back().id;

  // Loop over local cells
  double local_fission_rate = 0.0;
  for (auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& cell_matrices = unit_cell_matrices[cell.local_id];

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& sigma_f = xs.GetSigmaFission();

    // skip non-fissionable material
    if (not xs.IsFissionable())
      continue;

    // Loop over nodes
    const int num_nodes = transport_view.GetNumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.intV_shapeI(i);

      // Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
        local_fission_rate += sigma_f[g] * phi[uk_map + g] * IntV_ShapeI;
    } // for node
  } // for cell

  // Allreduce global production
  double global_fission_rate = 0.0;
  mpi_comm.all_reduce(local_fission_rate, global_fission_rate, mpi::op::sum<double>());

  return global_fission_rate;
}

void
ComputePrecursors(LBSProblem& lbs_problem)
{
  CALI_CXX_MARK_SCOPE("ComputePrecursors");

  const size_t J = lbs_problem.GetMaxPrecursorsPerMaterial();

  auto& precursor_new_local = lbs_problem.GetPrecursorsNewLocal();
  precursor_new_local.assign(precursor_new_local.size(), 0.0);

  auto& grid = lbs_problem.GetGrid();
  auto& groups = lbs_problem.GetGroups();
  auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto& phi_new_local = lbs_problem.GetPhiNewLocal();

  // Loop over cells
  for (const auto& cell : grid->local_cells)
  {
    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& transport_view = cell_transport_views[cell.local_id];
    const double cell_volume = transport_view.GetVolume();

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& precursors = xs.GetPrecursors();
    const auto& nu_delayed_sigma_f = xs.GetNuDelayedSigmaF();

    // Loop over precursors
    for (uint64_t j = 0; j < xs.GetNumPrecursors(); ++j)
    {
      size_t dof = cell.local_id * J + j;
      const auto& precursor = precursors[j];
      const double coeff = precursor.fractional_yield / precursor.decay_constant;

      // Loop over nodes
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        const size_t uk_map = transport_view.MapDOF(i, 0, 0);
        const double node_V_fraction = fe_values.intV_shapeI(i) / cell_volume;

        // Loop over groups
        for (unsigned int g = 0; g < groups.size(); ++g)
          precursor_new_local[dof] +=
            coeff * nu_delayed_sigma_f[g] * phi_new_local[uk_map + g] * node_V_fraction;
      } // for node i
    } // for precursor j
  } // for cell
}

} // namespace opensn
