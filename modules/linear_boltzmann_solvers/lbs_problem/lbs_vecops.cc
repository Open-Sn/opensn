// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/math/petsc_utils/petsc_utils.h"

namespace opensn
{

template <typename Functor>
int64_t
LBSVecOps::GroupsetScopedCopy(LBSProblem& lbs_problem,
                              unsigned int gsi,
                              unsigned int gss,
                              Functor func)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::GroupsetScopedCopy");

  const auto& grid = lbs_problem.GetGrid();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto num_moments = lbs_problem.GetNumMoments();

  int64_t idx = -1;
  for (const auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    for (std::size_t i = 0; i < cell.vertex_ids.size(); ++i)
    {
      for (unsigned int m = 0; m < num_moments; ++m)
      {
        auto mapped_idx = transport_view.MapDOF(i, m, gsi);
        for (unsigned int g = 0; g < gss; ++g)
        {
          ++idx;
          func(idx, mapped_idx + g);
        }
      }
    }
  }
  return idx;
}

void
LBSVecOps::SetPhiVectorScalarValues(LBSProblem& lbs_problem, PhiSTLOption phi_opt, double value)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::SetPhiVectorScalarValues");

  auto& phi = (phi_opt == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()
                                                 : lbs_problem.GetPhiOldLocal();
  const auto& grid = lbs_problem.GetGrid();
  const auto& sdm = lbs_problem.GetSpatialDiscretization();
  const auto& unknown_manager = lbs_problem.GetUnknownManager();

  const auto first_grp = 0;
  const long final_grp = static_cast<long>(lbs_problem.GetNumGroups() - 1);

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dof_map = static_cast<long>(sdm.MapDOFLocal(cell, i, unknown_manager, 0, 0));
      std::fill(phi.begin() + dof_map + first_grp, phi.begin() + dof_map + final_grp + 1, value);
    }
  }
}

void
LBSVecOps::ScalePhiVector(LBSProblem& lbs_problem, PhiSTLOption phi_opt, double value)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::ScalePhiVector");

  auto& phi = (phi_opt == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()
                                                 : lbs_problem.GetPhiOldLocal();
  auto& groupsets = lbs_problem.GetGroupsets();

  Scale(phi, value);
  for (auto& groupset : groupsets)
  {
    if (groupset.angle_agg)
    {
      if (phi_opt == PhiSTLOption::PHI_NEW)
      {
        auto psi = groupset.angle_agg->GetNewDelayedAngularDOFsAsSTLVector();
        Scale(psi, value);
        groupset.angle_agg->SetNewDelayedAngularDOFsFromSTLVector(psi);
      }
      else if (phi_opt == PhiSTLOption::PHI_OLD)
      {
        auto psi = groupset.angle_agg->GetOldDelayedAngularDOFsAsSTLVector();
        Scale(psi, value);
        groupset.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(psi);
      }
    }
  }
}

void
LBSVecOps::SetGSPETScVecFromPrimarySTLvector(LBSProblem& lbs_problem,
                                             const LBSGroupset& groupset,
                                             Vec dest,
                                             PhiSTLOption src)
{
  const auto& src_phi =
    (src == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  double* petsc_dest = nullptr;
  OpenSnPETScCall(VecGetArray(dest, &petsc_dest));
  int64_t index = GroupsetScopedCopy(lbs_problem,
                                     groupset.first_group,
                                     groupset.GetNumGroups(),
                                     [&](int64_t idx, size_t mapped_idx)
                                     { petsc_dest[idx] = src_phi[mapped_idx]; });
  if (groupset.angle_agg)
  {
    if (src == PhiSTLOption::PHI_NEW)
      groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, petsc_dest);
    else if (src == PhiSTLOption::PHI_OLD)
      groupset.angle_agg->AppendOldDelayedAngularDOFsToArray(index, petsc_dest);
  }
  OpenSnPETScCall(VecRestoreArray(dest, &petsc_dest));
}

void
LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(LBSProblem& lbs_problem,
                                             const LBSGroupset& groupset,
                                             Vec src,
                                             PhiSTLOption dest)
{
  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  const double* petsc_src = nullptr;
  OpenSnPETScCall(VecGetArrayRead(src, &petsc_src));
  int64_t index = GroupsetScopedCopy(lbs_problem,
                                     groupset.first_group,
                                     groupset.GetNumGroups(),
                                     [&](int64_t idx, size_t mapped_idx)
                                     { dest_phi[mapped_idx] = petsc_src[idx]; });
  if (groupset.angle_agg)
  {
    if (dest == PhiSTLOption::PHI_NEW)
      groupset.angle_agg->SetNewDelayedAngularDOFsFromArray(index, petsc_src);
    else if (dest == PhiSTLOption::PHI_OLD)
      groupset.angle_agg->SetOldDelayedAngularDOFsFromArray(index, petsc_src);
  }
  OpenSnPETScCall(VecRestoreArrayRead(src, &petsc_src));
}

void
LBSVecOps::SetGroupScopedPETScVecFromPrimarySTLvector(LBSProblem& lbs_problem,
                                                      unsigned int first_group_id,
                                                      unsigned int last_group_id,
                                                      Vec dest,
                                                      const std::vector<double>& src)
{
  double* petsc_dest = nullptr;
  OpenSnPETScCall(VecGetArray(dest, &petsc_dest));
  GroupsetScopedCopy(lbs_problem,
                     first_group_id,
                     last_group_id - first_group_id + 1,
                     [&](int64_t idx, size_t mapped_idx) { petsc_dest[idx] = src[mapped_idx]; });
  OpenSnPETScCall(VecRestoreArray(dest, &petsc_dest));
}

void
LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(LBSProblem& lbs_problem,
                                                      unsigned int first_group_id,
                                                      unsigned int last_group_id,
                                                      Vec src,
                                                      PhiSTLOption dest)
{
  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  const double* petsc_src = nullptr;
  OpenSnPETScCall(VecGetArrayRead(src, &petsc_src));
  GroupsetScopedCopy(lbs_problem,
                     first_group_id,
                     last_group_id - first_group_id + 1,
                     [&](int64_t idx, size_t mapped_idx)
                     { dest_phi[mapped_idx] = petsc_src[idx]; });
  OpenSnPETScCall(VecRestoreArrayRead(src, &petsc_src));
}

void
LBSVecOps::GSScopedCopyPrimarySTLvectors(LBSProblem& lbs_problem,
                                         const LBSGroupset& groupset,
                                         const std::vector<double>& src,
                                         std::vector<double>& dest)
{
  GroupsetScopedCopy(lbs_problem,
                     groupset.first_group,
                     groupset.GetNumGroups(),
                     [&](int64_t idx, size_t mapped_idx) { dest[mapped_idx] = src[mapped_idx]; });
}

void
LBSVecOps::GSScopedCopyPrimarySTLvectors(LBSProblem& lbs_problem,
                                         const LBSGroupset& groupset,
                                         PhiSTLOption src,
                                         PhiSTLOption dest)
{
  const auto& src_phi =
    (src == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  GroupsetScopedCopy(lbs_problem,
                     groupset.first_group,
                     groupset.GetNumGroups(),
                     [&](int64_t idx, size_t mapped_idx)
                     { dest_phi[mapped_idx] = src_phi[mapped_idx]; });
  if (groupset.angle_agg)
  {
    if (src == PhiSTLOption::PHI_NEW and dest == PhiSTLOption::PHI_OLD)
      groupset.angle_agg->SetDelayedPsiOld2New();
    else if (src == PhiSTLOption::PHI_OLD and dest == PhiSTLOption::PHI_NEW)
      groupset.angle_agg->SetDelayedPsiNew2Old();
  }
}

void
LBSVecOps::SetMultiGSPETScVecFromPrimarySTLvector(LBSProblem& lbs_problem,
                                                  const std::vector<unsigned int>& groupset_ids,
                                                  Vec x,
                                                  PhiSTLOption which_phi)
{
  auto& y = (which_phi == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()
                                                 : lbs_problem.GetPhiOldLocal();
  auto& groupsets = lbs_problem.GetGroupsets();
  double* x_ref = nullptr;
  OpenSnPETScCall(VecGetArray(x, &x_ref));

  for (auto gs_id : groupset_ids)
  {
    const auto& groupset = groupsets.at(gs_id);
    auto gsi = groupset.first_group;
    auto gsf = groupset.last_group;
    auto gss = gsf - gsi + 1;

    int64_t index = GroupsetScopedCopy(
      lbs_problem, gsi, gss, [&](int64_t idx, size_t mapped_idx) { x_ref[idx] = y[mapped_idx]; });
    if (groupset.angle_agg)
    {
      if (which_phi == PhiSTLOption::PHI_NEW)
        groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, x_ref);
      else if (which_phi == PhiSTLOption::PHI_OLD)
        groupset.angle_agg->AppendOldDelayedAngularDOFsToArray(index, x_ref);
    }
  }

  OpenSnPETScCall(VecRestoreArray(x, &x_ref));
}

void
LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec(LBSProblem& lbs_problem,
                                                  const std::vector<unsigned int>& groupset_ids,
                                                  Vec x,
                                                  PhiSTLOption which_phi)
{
  auto& y = (which_phi == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()
                                                 : lbs_problem.GetPhiOldLocal();
  auto& groupsets = lbs_problem.GetGroupsets();
  const double* x_ref = nullptr;
  OpenSnPETScCall(VecGetArrayRead(x, &x_ref));

  for (auto gs_id : groupset_ids)
  {
    const auto& groupset = groupsets.at(gs_id);
    auto gsi = groupset.first_group;
    auto gsf = groupset.last_group;
    auto gss = gsf - gsi + 1;

    int64_t index = GroupsetScopedCopy(
      lbs_problem, gsi, gss, [&](int64_t idx, size_t mapped_idx) { y[mapped_idx] = x_ref[idx]; });
    if (groupset.angle_agg)
    {
      if (which_phi == PhiSTLOption::PHI_NEW)
        groupset.angle_agg->SetNewDelayedAngularDOFsFromArray(index, x_ref);
      else if (which_phi == PhiSTLOption::PHI_OLD)
        groupset.angle_agg->SetOldDelayedAngularDOFsFromArray(index, x_ref);
    }
  }

  OpenSnPETScCall(VecRestoreArrayRead(x, &x_ref));
}

} // namespace opensn
