// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/spatial_weight_function.h"
#include "framework/runtime.h"
#include <cassert>

namespace opensn
{

bool
UseDelayedNeutronProduction(const LBSProblem& lbs_problem)
{
  return lbs_problem.GetNumPrecursors() > 0 and lbs_problem.GetOptions().use_precursors;
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
  const bool adjoint = lbs_problem.GetOptions().adjoint;

  const auto first_grp = 0;
  const auto last_grp = lbs_problem.GetNumGroups() - 1;

  // Delayed production of a cell summed over the destination group, per unit flux in group g'
  std::vector<double> delayed_production(last_grp + 1, 0.0);

  // Loop over local cells
  double local_production = 0.0;
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& cell_matrices = unit_cell_matrices[cell.local_id];

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& F = xs.GetProductionMatrix(adjoint);

    if (not xs.IsFissionable())
      continue;

    // When delayed-neutron data is included, the production matrix is prompt-only and delayed
    // production is added separately. If the production matrix changes to include delayed
    // production, adjust this sum to prevent double counting. In adjoint mode F is transposed, and
    // so is the delayed production.
    if (include_delayed_production)
      for (unsigned int gp = first_grp; gp <= last_grp; ++gp)
      {
        delayed_production[gp] = 0.0;
        for (unsigned int g = first_grp; g <= last_grp; ++g)
          delayed_production[gp] += adjoint ? ComputeDelayedFissionProduction(xs, gp, g)
                                            : ComputeDelayedFissionProduction(xs, g, gp);
      }

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
      }

      if (include_delayed_production)
        for (unsigned int gp = first_grp; gp <= last_grp; ++gp)
          local_production += delayed_production[gp] * phi[uk_map + gp] * IntV_ShapeI;
    } // for node
  } // for cell

  // Report a physical rate: curvilinear unit integrals omit the angular extent (2*pi in RZ).
  local_production *= IntegralMeasureScale(grid->GetCoordinateSystem());

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
  for (const auto& cell : grid->GetLocalCells())
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

  // Report a physical rate: curvilinear unit integrals omit the angular extent (2*pi in RZ).
  local_fission_rate *= IntegralMeasureScale(grid->GetCoordinateSystem());

  // Allreduce global fission rate
  double global_fission_rate = 0.0;
  mpi_comm.all_reduce(local_fission_rate, global_fission_rate, mpi::op::sum<double>());

  return global_fission_rate;
}

void
ComputePrecursors(LBSProblem& lbs_problem)
{
  // Equilibrium precursor concentrations C_j = gamma_j / lambda_j * sum_g nu_d sigma_f,g phi_g,
  // stored per spatial node like the scalar flux (index node * J + j).
  const auto J = lbs_problem.GetMaxPrecursorsPerMaterial();

  auto& precursor_new_local = lbs_problem.GetPrecursorsNewLocal();
  precursor_new_local.assign(precursor_new_local.size(), 0.0);

  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  const auto& phi_new_local = lbs_problem.GetPhiNewLocal();

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& xs = transport_view.GetXS();
    const auto& precursors = xs.GetPrecursors();
    if (precursors.empty())
      continue;
    const auto& nu_delayed_sigma_f = xs.GetNuDelayedSigmaF();

    for (int i = 0; i < transport_view.GetNumNodes(); ++i)
    {
      const auto uk_map = transport_view.MapDOF(i, 0, 0);
      double delayed_production = 0.0;
      for (unsigned int g = 0; g < lbs_problem.GetNumGroups(); ++g)
        delayed_production += nu_delayed_sigma_f[g] * phi_new_local[uk_map + g];

      const size_t node_base = discretization.MapDOFLocal(cell, i) * J;
      for (unsigned int j = 0; j < precursors.size(); ++j)
      {
        const auto& precursor = precursors[j];
        assert(precursor.decay_constant > 0.0 &&
               "ComputePrecursors encountered non-positive precursor decay constant.");
        precursor_new_local[node_base + j] =
          precursor.fractional_yield / precursor.decay_constant * delayed_production;
      }
    } // for node i
  } // for cell
}

} // namespace opensn
