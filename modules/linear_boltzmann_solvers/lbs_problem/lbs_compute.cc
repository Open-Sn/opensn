// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "caliper/cali.h"

namespace opensn
{

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
