// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules//linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/math/linear_solver/linear_matrix_action_Ax.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <petscksp.h>
#include <iomanip>

namespace opensn
{

void
AGSLinearSolver::Solve()
{
  CALI_CXX_MARK_SCOPE("AGSLinearSolver::Solve");

  auto ags_context_ptr = std::dynamic_pointer_cast<AGSContext>(context_ptr_);
  auto grid_ptr = GetCurrentMesh();
  auto& lbs_solver = ags_context_ptr->lbs_solver_;
  auto& wgs_solvers = ags_context_ptr->wgs_solvers_;
  auto& cell_transport_views = lbs_solver.GetCellTransportViews();
  auto& groupsets = lbs_solver.Groupsets();
  auto& phi_new = lbs_solver.PhiNewLocal();
  auto num_moments = lbs_solver.NumMoments();
  std::vector<double> phi_old(lbs_solver.PhiOldLocal().size(), 0.0);

  auto piecewise_change =
    [&grid_ptr, &cell_transport_views, &groupsets, num_moments, &phi_new, &phi_old]() {
      double pw_change = 0.0;
      for (const auto& cell : grid_ptr->local_cells)
      {
        auto& transport_view = cell_transport_views[cell.local_id_];
        for (int i = 0; i < cell.vertex_ids_.size(); ++i)
        {
          for (auto groupset : groupsets)
          {
            int gsi = groupset.groups_.front().id_;
            for (int g = 0; g < groupset.groups_.size(); ++g)
            {
              size_t m0g_idx = transport_view.MapDOF(i, 0, gsi + g);
              double max_phi = std::max(fabs(phi_new[m0g_idx]), fabs(phi_old[m0g_idx]));
              for (int m = 0; m < num_moments; ++m)
              {
                size_t mng_idx = transport_view.MapDOF(i, m, gsi + g);
                double delta_phi = std::fabs(phi_new[mng_idx] - phi_old[mng_idx]);
                if (max_phi >= std::numeric_limits<double>::min())
                  pw_change = std::max(delta_phi / max_phi, pw_change);
                else
                  pw_change = std::max(delta_phi, pw_change);
              } // for g
            }   // for m
          }     // for groupset
        }       // for i
      }         // for c
      double global_pw_change = 0.0;
      mpi_comm.all_reduce<double>(pw_change, global_pw_change, mpi::op::max<double>());
      return global_pw_change;
    };

  // Save qmoms to be restored after each iteration. This is necessary for multiple ags iterations
  // to function and for keigen-value problems
  const auto saved_qmoms = lbs_solver.QMomentsLocal();

  double pw_change_prev = 1.0;
  bool converged = false;
  for (int iter = 0; iter < maximum_iterations_; ++iter)
  {
    for (auto& solver : wgs_solvers)
    {
      solver->Setup();
      solver->Solve();
    }

    double pw_change = piecewise_change();
    double rho = (iter == 0) ? 0.0 : sqrt(pw_change / pw_change_prev);
    pw_change_prev = pw_change;

    if (pw_change < std::max(tolerance_ * (1.0 - rho), 1.0e-10))
      converged = true;
    else
      phi_old = lbs_solver.PhiNewLocal();

    lbs_solver.QMomentsLocal() = saved_qmoms; // Restore qmoms

    // Write restart data
    if (lbs_solver.RestartsEnabled() and lbs_solver.TriggerRestartDump() and
        lbs_solver.Options().enable_ags_restart_write)
    {
      lbs_solver.WriteRestartData();
    }

    if (verbose_)
    {
      std::stringstream iter_stats;
      iter_stats << " AGS:"
                 << " Iteration = " << std::left << std::setw(5) << iter
                 << " Point-wise change = " << std::left << std::setw(14) << pw_change
                 << " Spectral-radius estimate = " << std::left << std::setw(10) << rho;
      if (converged)
        iter_stats << " CONVERGED";
      log.Log() << iter_stats.str();
    }

    if (converged)
      break;
  }

  // If restarts are enabled, always write a restart dump upon convergence or when we reach the
  // iteration limit
  if (lbs_solver.RestartsEnabled() && lbs_solver.Options().enable_ags_restart_write)
    lbs_solver.WriteRestartData();
}

} // namespace opensn
