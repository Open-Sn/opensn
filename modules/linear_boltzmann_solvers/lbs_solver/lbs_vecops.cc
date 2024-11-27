// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_vecops.h"

namespace opensn
{

void
LBSVecOps::SetPhiVectorScalarValues(LBSSolver& lbs_solver, PhiSTLOption phi_opt, double value)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::SetPhiVectorScalarValues");

  auto& phi =
    (phi_opt == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  auto& grid = lbs_solver.Grid();
  auto& groups = lbs_solver.Groups();
  auto& sdm = lbs_solver.SpatialDiscretization();
  auto& unknown_manager = lbs_solver.UnknownManager();

  const size_t first_grp = groups.front().id;
  const size_t final_grp = groups.back().id;

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dof_map = sdm.MapDOFLocal(cell, i, unknown_manager, 0, 0);
      std::fill(phi.begin() + dof_map + first_grp, phi.begin() + dof_map + final_grp + 1, value);
    }
  }
}

void
LBSVecOps::ScalePhiVector(LBSSolver& lbs_solver, PhiSTLOption phi_opt, double value)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::ScalePhiVector");

  auto& phi =
    (phi_opt == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  auto& groupsets = lbs_solver.Groupsets();

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
LBSVecOps::SetGSPETScVecFromPrimarySTLvector(LBSSolver& lbs_solver,
                                             const LBSGroupset& groupset,
                                             Vec dest,
                                             PhiSTLOption src)
{
  const auto& src_phi =
    (src == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  double* petsc_dest;
  VecGetArray(dest, &petsc_dest);
  int64_t index = GroupsetScopedCopy(lbs_solver,
                                     groupset.groups.front().id,
                                     groupset.groups.size(),
                                     [&](int64_t idx, size_t mapped_idx)
                                     { petsc_dest[idx] = src_phi[mapped_idx]; });
  if (groupset.angle_agg)
  {
    if (src == PhiSTLOption::PHI_NEW)
      groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, petsc_dest);
    else if (src == PhiSTLOption::PHI_OLD)
      groupset.angle_agg->AppendOldDelayedAngularDOFsToArray(index, petsc_dest);
  }
  VecRestoreArray(dest, &petsc_dest);
}

void
LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(LBSSolver& lbs_solver,
                                             const LBSGroupset& groupset,
                                             Vec src,
                                             PhiSTLOption dest)
{
  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  const double* petsc_src;
  VecGetArrayRead(src, &petsc_src);
  int64_t index = GroupsetScopedCopy(lbs_solver,
                                     groupset.groups.front().id,
                                     groupset.groups.size(),
                                     [&](int64_t idx, size_t mapped_idx)
                                     { dest_phi[mapped_idx] = petsc_src[idx]; });
  if (groupset.angle_agg)
  {
    if (dest == PhiSTLOption::PHI_NEW)
      groupset.angle_agg->SetNewDelayedAngularDOFsFromArray(index, petsc_src);
    else if (dest == PhiSTLOption::PHI_OLD)
      groupset.angle_agg->SetOldDelayedAngularDOFsFromArray(index, petsc_src);
  }
  VecRestoreArrayRead(src, &petsc_src);
}

void
LBSVecOps::SetGroupScopedPETScVecFromPrimarySTLvector(LBSSolver& lbs_solver,
                                                      int first_group_id,
                                                      int last_group_id,
                                                      Vec dest,
                                                      const std::vector<double>& src)
{
  double* petsc_dest;
  VecGetArray(dest, &petsc_dest);
  GroupsetScopedCopy(lbs_solver,
                     first_group_id,
                     last_group_id - first_group_id + 1,
                     [&](int64_t idx, size_t mapped_idx) { petsc_dest[idx] = src[mapped_idx]; });
  VecRestoreArray(dest, &petsc_dest);
}

void
LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(
  LBSSolver& lbs_solver, int first_group_id, int last_group_id, Vec src, PhiSTLOption dest)
{
  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  const double* petsc_src;
  VecGetArrayRead(src, &petsc_src);
  GroupsetScopedCopy(lbs_solver,
                     first_group_id,
                     last_group_id - first_group_id + 1,
                     [&](int64_t idx, size_t mapped_idx)
                     { dest_phi[mapped_idx] = petsc_src[idx]; });
  VecRestoreArrayRead(src, &petsc_src);
}

void
LBSVecOps::GSScopedCopyPrimarySTLvectors(LBSSolver& lbs_solver,
                                         const LBSGroupset& groupset,
                                         const std::vector<double>& src,
                                         std::vector<double>& dest)
{
  GroupsetScopedCopy(lbs_solver,
                     groupset.groups.front().id,
                     groupset.groups.size(),
                     [&](int64_t idx, size_t mapped_idx) { dest[mapped_idx] = src[mapped_idx]; });
}

void
LBSVecOps::GSScopedCopyPrimarySTLvectors(LBSSolver& lbs_solver,
                                         const LBSGroupset& groupset,
                                         PhiSTLOption src,
                                         PhiSTLOption dest)
{
  const auto& src_phi =
    (src == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  int64_t index = GroupsetScopedCopy(lbs_solver,
                                     groupset.groups.front().id,
                                     groupset.groups.size(),
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
LBSVecOps::SetMultiGSPETScVecFromPrimarySTLvector(LBSSolver& lbs_solver,
                                                  const std::vector<int>& groupset_ids,
                                                  Vec x,
                                                  PhiSTLOption which_phi)
{
  auto& y =
    (which_phi == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  auto& groupsets = lbs_solver.Groupsets();
  double* x_ref;
  VecGetArray(x, &x_ref);

  for (int gs_id : groupset_ids)
  {
    const auto& groupset = groupsets.at(gs_id);
    int gsi = groupset.groups.front().id;
    int gsf = groupset.groups.back().id;
    int gss = gsf - gsi + 1;

    int64_t index = GroupsetScopedCopy(
      lbs_solver, gsi, gss, [&](int64_t idx, size_t mapped_idx) { x_ref[idx] = y[mapped_idx]; });
    if (groupset.angle_agg)
    {
      if (which_phi == PhiSTLOption::PHI_NEW)
        groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, x_ref);
      else if (which_phi == PhiSTLOption::PHI_OLD)
        groupset.angle_agg->AppendOldDelayedAngularDOFsToArray(index, x_ref);
    }
  }

  VecRestoreArray(x, &x_ref);
}

void
LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec(LBSSolver& lbs_solver,
                                                  const std::vector<int>& groupset_ids,
                                                  Vec x,
                                                  PhiSTLOption which_phi)
{
  auto& y =
    (which_phi == PhiSTLOption::PHI_NEW) ? lbs_solver.PhiNewLocal() : lbs_solver.PhiOldLocal();
  auto& groupsets = lbs_solver.Groupsets();
  const double* x_ref;
  VecGetArrayRead(x, &x_ref);

  for (int gs_id : groupset_ids)
  {
    const auto& groupset = groupsets.at(gs_id);
    int gsi = groupset.groups.front().id;
    int gsf = groupset.groups.back().id;
    int gss = gsf - gsi + 1;

    int64_t index = GroupsetScopedCopy(
      lbs_solver, gsi, gss, [&](int64_t idx, size_t mapped_idx) { y[mapped_idx] = x_ref[idx]; });
    if (groupset.angle_agg)
    {
      if (which_phi == PhiSTLOption::PHI_NEW)
        groupset.angle_agg->SetNewDelayedAngularDOFsFromArray(index, x_ref);
      else if (which_phi == PhiSTLOption::PHI_OLD)
        groupset.angle_agg->SetOldDelayedAngularDOFsFromArray(index, x_ref);
    }
  }

  VecRestoreArrayRead(x, &x_ref);
}

} // namespace opensn
