// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules//linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_solver.h"
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
AGSSolver::Solve()
{
  CALI_CXX_MARK_SCOPE("AGSSolver::Solve");

  std::fill(phi_old_.begin(), phi_old_.end(), 0.0);

  // Save qmoms to be restored after each iteration. This is necessary for multiple ags iterations
  // to function and for keigen-value problems
  const auto saved_qmoms = lbs_solver_.QMomentsLocal();

  double pw_change_prev = 1.0;
  bool converged = false;
  for (int iter = 0; iter < max_iterations_; ++iter)
  {
    for (auto& solver : wgs_solvers_)
    {
      solver->Setup();
      solver->Solve();
    }

    std::stringstream iter_stats;
    iter_stats << program_timer.GetTimeString() << " AGS Iteration ";

    if (lbs_solver_.Options().ags_pointwise_convergence)
    {
      double pw_change = ComputePointwisePhiChange();
      double rho = (iter == 0) ? 0.0 : sqrt(pw_change / pw_change_prev);
      pw_change_prev = pw_change;

      iter_stats << std::left << std::setw(5) << iter << " Point-wise change " << std::left
                 << std::setw(14) << pw_change << " Spectral-radius estimate " << std::left
                 << std::setw(10) << rho;

      if (pw_change < std::max(tolerance_ * (1.0 - rho), 1.0e-10))
      {
        converged = true;
        iter_stats << " CONVERGED";
      }
    }
    else
    {
      double norm = ComputeL2PhiChange();

      iter_stats << std::left << std::setw(5) << iter << " Error Norm " << std::left
                 << std::setw(14) << norm;

      if (norm < tolerance_)
      {
        converged = true;
        iter_stats << " CONVERGED";
      }
    }

    if (verbose_)
      log.Log() << iter_stats.str();

    // Restore qmoms
    lbs_solver_.QMomentsLocal() = saved_qmoms;

    // Write restart data
    if (lbs_solver_.RestartsEnabled() and lbs_solver_.TriggerRestartDump() and
        lbs_solver_.Options().enable_ags_restart_write)
    {
      lbs_solver_.WriteRestartData();
    }

    if (converged)
      break;
    else
      phi_old_ = lbs_solver_.PhiNewLocal();
  }

  // If restarts are enabled, always write a restart dump upon convergence or when we reach the
  // iteration limit
  if (lbs_solver_.RestartsEnabled() && lbs_solver_.Options().enable_ags_restart_write)
    lbs_solver_.WriteRestartData();
}

double
AGSSolver::ComputeL2PhiChange() const
{
  double norm = 0.0;
  auto& phi_new = lbs_solver_.PhiNewLocal();
  for (int i = 0; i < phi_new.size(); ++i)
  {
    double val = phi_new[i] - phi_old_[i];
    norm += val * val;
  }

  double global_norm = 0.0;
  mpi_comm.all_reduce<double>(norm, global_norm, mpi::op::sum<double>());
  global_norm = std::sqrt(global_norm);

  return global_norm;
}

double
AGSSolver::ComputePointwisePhiChange() const
{
  auto grid_ptr = GetCurrentMesh();
  auto& cell_transport_views = lbs_solver_.GetCellTransportViews();
  auto& groupsets = lbs_solver_.Groupsets();
  auto& phi_new = lbs_solver_.PhiNewLocal();
  auto num_moments = lbs_solver_.NumMoments();

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
          double max_phi = std::max(fabs(phi_new[m0g_idx]), fabs(phi_old_[m0g_idx]));
          for (int m = 0; m < num_moments; ++m)
          {
            size_t mng_idx = transport_view.MapDOF(i, m, gsi + g);
            double delta_phi = std::fabs(phi_new[mng_idx] - phi_old_[mng_idx]);
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
}

} // namespace opensn