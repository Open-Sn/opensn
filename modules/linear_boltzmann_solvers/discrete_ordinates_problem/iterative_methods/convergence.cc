// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/convergence.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/utils/error.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh/mesh.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace opensn
{

double
ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  int groupset_id,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old)
{
  return ComputePointwisePhiChange(lbs_problem, std::vector<int>{groupset_id}, opt_phi_old);
}

double
ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old)
{
  return ComputePointwisePhiChange(lbs_problem, std::vector<int>{}, opt_phi_old);
}

double
ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  std::vector<int> groupset_ids,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old)
{
  if (groupset_ids.empty())
  {
    groupset_ids.resize(lbs_problem.GetNumGroupsets());
    std::iota(groupset_ids.begin(), groupset_ids.end(), 0);
  }

  const auto& phi_new = lbs_problem.GetPhiNewLocal();
  const std::vector<double>& phi_old =
    opt_phi_old.has_value() ? opt_phi_old.value().get() : lbs_problem.GetPhiOldLocal();

  auto grid_ptr = lbs_problem.GetGrid();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto num_moments = lbs_problem.GetNumMoments();

  double pw_change = 0.0;
  for (const auto& cell : grid_ptr->GetLocalCells())
  {
    const auto& transport_view = cell_transport_views[cell->local_id];
    for (auto i = 0; i < cell->vertex_ids.size(); ++i)
    {
      for (auto id : groupset_ids)
      {
        const auto& groupset = lbs_problem.GetGroupset(id);
        auto gsi = groupset.first_group;
        for (unsigned int g = 0; g < groupset.GetNumGroups(); ++g)
        {
          // A higher moment that is negligible next to the scalar flux carries only round-off, so
          // its change relative to itself never settles. Each moment's change is therefore taken
          // relative to the larger of its own magnitude and the scalar flux at the same node.
          const auto phi0_idx = transport_view.MapDOF(i, 0, gsi + g);
          const double phi0_scale =
            std::max(std::fabs(phi_new[phi0_idx]), std::fabs(phi_old[phi0_idx]));
          for (unsigned int m = 0; m < num_moments; ++m)
          {
            auto mng_idx = transport_view.MapDOF(i, m, gsi + g);
            double delta_phi = std::fabs(phi_new[mng_idx] - phi_old[mng_idx]);
            const double scale =
              std::max({phi0_scale, std::fabs(phi_new[mng_idx]), std::fabs(phi_old[mng_idx])});
            if (scale >= std::numeric_limits<double>::min())
              pw_change = std::max(delta_phi / scale, pw_change);
            else
              pw_change = std::max(delta_phi, pw_change);
          } // for m
        } // for g
      } // for id
    } // for i
  } // for cell

  double global_pw_change = 0.0;
  mpi_comm.all_reduce<double>(pw_change, global_pw_change, mpi::op::max<double>());

  return global_pw_change;
}

double
ComputePointwisePsiChange(DiscreteOrdinatesProblem& do_problem,
                          const LBSGroupset& groupset,
                          const std::vector<double>& psi_new,
                          const std::vector<double>& psi_old)
{
  OpenSnLogicalErrorIf(psi_new.size() != psi_old.size(),
                       "ComputePointwisePsiChange: new and old delayed angular unknowns differ in "
                       "size.");

  double pw_change = 0.0;
  if (not psi_new.empty())
  {
    const auto& grid = *do_problem.GetGrid();
    const auto& sdm = do_problem.GetSpatialDiscretization();
    const auto& transport_views = do_problem.GetCellTransportViews();
    const auto& phi_new = do_problem.GetPhiNewLocal();
    const auto& phi_old = do_problem.GetPhiOldLocal();
    const auto first_group = groupset.first_group;

    // A lagged angular flux that is negligible next to the scalar flux where it is read, such as a
    // grazing direction behind a cycle, can hold only round-off, so its change relative to itself
    // never settles. Each entry's change is therefore taken relative to the larger of its own
    // magnitude and that scalar flux, which is the isotropic angular flux because the quadrature
    // weights sum to one.
    const auto phi_scale = [&](std::uint32_t cell_local_id,
                               unsigned int face,
                               unsigned int face_node,
                               unsigned int group_idx)
    {
      const auto& cell = *grid.GetLocalCells()[cell_local_id];
      const auto node = sdm.GetCellMapping(cell).MapFaceNode(face, face_node);
      const auto dof = transport_views[cell_local_id].MapDOF(node, 0, first_group + group_idx);
      return std::max(std::fabs(phi_new[dof]), std::fabs(phi_old[dof]));
    };
    OpenSnLogicalErrorIf(not groupset.angle_agg,
                         "ComputePointwisePsiChange: groupset has no angle aggregation.");
    const auto scales = groupset.angle_agg->MapNodalValuesToDelayedAngularDOFs(phi_scale);
    OpenSnLogicalErrorIf(scales.size() != psi_new.size(),
                         "ComputePointwisePsiChange: delayed angular unknowns differ in size from "
                         "the groupset's delayed storage.");

    for (size_t i = 0; i < psi_new.size(); ++i)
    {
      const double delta = std::fabs(psi_new[i] - psi_old[i]);
      const double scale = std::max({scales[i], std::fabs(psi_new[i]), std::fabs(psi_old[i])});
      if (scale >= std::numeric_limits<double>::min())
        pw_change = std::max(delta / scale, pw_change);
      else
        pw_change = std::max(delta, pw_change);
    }
  }

  double global_pw_change = 0.0;
  mpi_comm.all_reduce<double>(pw_change, global_pw_change, mpi::op::max<double>());

  return global_pw_change;
}

double
ComputeL2PhiChange(LBSProblem& lbs_problem,
                   std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old)
{
  double norm = 0.0;
  const auto& phi_new = lbs_problem.GetPhiNewLocal();
  const std::vector<double>& phi_old =
    opt_phi_old.has_value() ? opt_phi_old.value().get() : lbs_problem.GetPhiOldLocal();

  for (auto i = 0; i < phi_new.size(); ++i)
  {
    double val = phi_new[i] - phi_old[i];
    norm += val * val;
  }

  double global_norm = 0.0;
  mpi_comm.all_reduce<double>(norm, global_norm, mpi::op::sum<double>());

  return std::sqrt(global_norm);
}

double
ConvergenceRateEstimator::Update(const double change)
{
  ++num_changes_;
  expected_change_ = rate_ * previous_change_;
  if (num_changes_ == 1)
    rate_ = 0.0;
  else
  {
    const double ratio = previous_change_ > 0.0 ? change / previous_change_ : 1.0;
    rate_ = num_changes_ == 2 ? ratio : std::max(ratio, previous_ratio_);
    rate_ = std::clamp(rate_, 0.0, 1.0);
    previous_ratio_ = ratio;
  }
  previous_change_ = change;
  return rate_;
}

bool
ConvergenceRateEstimator::Passes(const double change,
                                 const double tolerance,
                                 const double roundoff_change) const
{
  return change < roundoff_change or change < tolerance * (1.0 - rate_);
}

bool
ConvergenceRateEstimator::Converged(const double tolerance, const double roundoff_change) const
{
  return Passes(std::max(previous_change_, expected_change_), tolerance, roundoff_change);
}

} // namespace opensn
