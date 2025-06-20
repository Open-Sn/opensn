// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/convergence.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <numeric>

namespace opensn
{

double
ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  int groupset_id,
  std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_phi_old)
{
  return ComputePointwisePhiChange(lbs_problem, std::vector<int>{groupset_id}, opt_phi_old);
}

double
ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_phi_old)
{
  return ComputePointwisePhiChange(lbs_problem, std::vector<int>{}, opt_phi_old);
}

double
ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  std::vector<int> groupset_ids,
  std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_phi_old)
{
  if (groupset_ids.empty())
  {
    groupset_ids.resize(lbs_problem.GetGroupsets().size());
    std::iota(groupset_ids.begin(), groupset_ids.end(), 0);
  }

  auto& phi_new = lbs_problem.GetPhiNewLocal();
  auto& phi_old =
    opt_phi_old.has_value() ? opt_phi_old.value().get() : lbs_problem.GetPhiOldLocal();

  auto grid_ptr = lbs_problem.GetGrid();
  auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto num_moments = lbs_problem.GetNumMoments();

  double pw_change = 0.0;
  for (auto gsi : groupset_ids)
  {
    auto& groupset = lbs_problem.GetGroupsets()[gsi];
    for (const auto& cell : grid_ptr->local_cells)
    {
      auto& transport_view = cell_transport_views[gsi][cell.local_id];
      for (auto i = 0; i < cell.vertex_ids.size(); ++i)
      {
        for (auto g = 0; g < groupset.groups.size(); ++g)
        {
          auto m0g_idx = transport_view.MapDOF(i, 0, g);
          double max_phi = std::max(fabs(phi_new[gsi][m0g_idx]), fabs(phi_old[gsi][m0g_idx]));
          for (auto m = 0; m < num_moments; ++m)
          {
            auto mng_idx = transport_view.MapDOF(i, m, g);
            double delta_phi = std::fabs(phi_new[gsi][mng_idx] - phi_old[gsi][mng_idx]);
            if (max_phi >= std::numeric_limits<double>::min())
              pw_change = std::max(delta_phi / max_phi, pw_change);
            else
              pw_change = std::max(delta_phi, pw_change);
          } // for m
        }   // for g
      }     // for i
    }       // for cell
  }         // for id

  double global_pw_change = 0.0;
  mpi_comm.all_reduce<double>(pw_change, global_pw_change, mpi::op::max<double>());

  return global_pw_change;
}

double
ComputeL2PhiChange(
  LBSProblem& lbs_problem,
  std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_phi_old)
{
  double norm = 0.0;

  auto& phi_new = lbs_problem.GetPhiNewLocal();
  auto& phi_old =
    opt_phi_old.has_value() ? opt_phi_old.value().get() : lbs_problem.GetPhiOldLocal();

  for (std::size_t gsi = 0; gsi < lbs_problem.GetGroupsets().size(); ++gsi)
  {
    for (std::size_t i = 0; i < phi_new[gsi].size(); ++i)
    {
      double val = phi_new[gsi][i] - phi_old[gsi][i];
      norm += val * val;
    }
  }

  double global_norm = 0.0;
  mpi_comm.all_reduce<double>(norm, global_norm, mpi::op::sum<double>());

  return std::sqrt(global_norm);
}

} // namespace opensn
