// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include <cassert>

namespace opensn
{

bool
UseDelayedNeutronProduction(const LBSProblem& lbs_problem)
{
  return lbs_problem.GetNumPrecursors() > 0 and
         (not lbs_problem.IsTimeDependent() or lbs_problem.GetOptions().use_precursors);
}

double
ComputeDelayedFissionProduction(const MultiGroupXS& xs,
                                const unsigned int to_group,
                                const unsigned int from_group)
{
  if (xs.GetNumPrecursors() == 0)
    return 0.0;

  const auto& nu_delayed_sigma_f = xs.GetNuDelayedSigmaF();
  if (nu_delayed_sigma_f.empty())
    return 0.0;

  double delayed_production = 0.0;
  for (const auto& precursor : xs.GetPrecursors())
    delayed_production += precursor.emission_spectrum[to_group] * precursor.fractional_yield *
                          nu_delayed_sigma_f[from_group];

  return delayed_production;
}

double
ComputeFissionProduction(const LBSProblem& lbs_problem, const std::vector<double>& phi)
{

  const auto& grid = lbs_problem.GetGrid();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  assert(phi.size() == lbs_problem.GetPhiNewLocal().size() &&
         "ComputeFissionProduction size mismatch.");
  const bool include_delayed_production = UseDelayedNeutronProduction(lbs_problem);

  const auto first_grp = 0;
  const auto last_grp = lbs_problem.GetNumGroups() - 1;

  // Loop over local cells
  double local_production = 0.0;
  for (auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& cell_matrices = unit_cell_matrices[cell.local_id];

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& F = xs.GetProductionMatrix();

    if (not xs.IsFissionable())
      continue;

    // Loop over nodes
    const int num_nodes = transport_view.GetNumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const auto uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.intV_shapeI(i);

      // Loop over groups
      for (unsigned int g = first_grp; g <= last_grp; ++g)
      {
        const auto& prod = F[g];
        for (unsigned int gp = 0; gp <= last_grp; ++gp)
          local_production += prod[gp] * phi[uk_map + gp] * IntV_ShapeI;

        // When delayed-neutron data is included, the production matrix is prompt-only
        // and delayed production is added separately. If the production matrix changes
        // to include delayed production, adjust this sum to prevent double counting.
        if (include_delayed_production)
          for (unsigned int gp = 0; gp <= last_grp; ++gp)
            local_production +=
              ComputeDelayedFissionProduction(xs, g, gp) * phi[uk_map + gp] * IntV_ShapeI;
      }
    } // for node
  } // for cell

  // Allreduce global production
  double global_production = 0.0;
  mpi_comm.all_reduce(local_production, global_production, mpi::op::sum<double>());

  return global_production;
}

double
ComputeFissionRate(const LBSProblem& lbs_problem, const std::vector<double>& phi)
{

  const auto& grid = lbs_problem.GetGrid();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  assert(phi.size() == lbs_problem.GetPhiNewLocal().size() && "ComputeFissionRate size mismatch.");

  const auto first_grp = 0;
  const auto last_grp = lbs_problem.GetNumGroups() - 1;

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
    for (auto i = 0; i < num_nodes; ++i)
    {
      const auto uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.intV_shapeI(i);

      // Loop over groups
      for (unsigned int g = first_grp; g <= last_grp; ++g)
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

  const auto J = lbs_problem.GetMaxPrecursorsPerMaterial();

  auto& precursor_new_local = lbs_problem.GetPrecursorsNewLocal();
  precursor_new_local.assign(precursor_new_local.size(), 0.0);

  const auto& grid = lbs_problem.GetGrid();
  const auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  const auto& phi_new_local = lbs_problem.GetPhiNewLocal();

  // Loop over cells
  for (const auto& cell : grid->local_cells)
  {
    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& transport_view = cell_transport_views[cell.local_id];
    const double cell_volume = transport_view.GetVolume();
    assert(cell_volume > 0.0 && "ComputePrecursors encountered non-positive cell volume.");

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& precursors = xs.GetPrecursors();
    const auto& nu_delayed_sigma_f = xs.GetNuDelayedSigmaF();

    // Loop over precursors
    for (unsigned int j = 0; j < xs.GetNumPrecursors(); ++j)
    {
      size_t dof = cell.local_id * J + j;
      const auto& precursor = precursors[j];
      assert(precursor.decay_constant > 0.0 &&
             "ComputePrecursors encountered non-positive precursor decay constant.");
      const double coeff = precursor.fractional_yield / precursor.decay_constant;

      // Loop over nodes
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        const auto uk_map = transport_view.MapDOF(i, 0, 0);
        const double node_V_fraction = fe_values.intV_shapeI(i) / cell_volume;

        // Loop over groups
        for (unsigned int g = 0; g < lbs_problem.GetNumGroups(); ++g)
          precursor_new_local[dof] +=
            coeff * nu_delayed_sigma_f[g] * phi_new_local[uk_map + g] * node_V_fraction;
      } // for node i
    } // for precursor j
  } // for cell
}

} // namespace opensn
