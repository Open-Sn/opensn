// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/math/petsc_utils/petsc_utils.h"

namespace opensn
{

namespace
{

/**
 * Calls `func(idx, mapped_idx)` for every local node, moment, and group in
 * [first_group, first_group + num_groups), with `idx` counting entries in [node][moment][group]
 * order. Unlike LBSVecOps::GroupsetScopedCopy, the groups may belong to several groupsets, so each
 * group is mapped separately.
 */
template <typename Functor>
void
GroupRangeCopy(const LBSProblem& lbs_problem,
               unsigned int first_group,
               unsigned int num_groups,
               Functor func)
{
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  const auto num_moments = lbs_problem.GetNumMoments();

  int64_t idx = -1;
  for (const auto& cell : lbs_problem.GetGrid()->GetLocalCells())
  {
    const auto& transport_view = cell_transport_views[cell->local_id];
    const auto num_nodes = static_cast<std::size_t>(transport_view.GetNumNodes());
    for (std::size_t i = 0; i < num_nodes; ++i)
      for (unsigned int m = 0; m < num_moments; ++m)
        for (unsigned int g = 0; g < num_groups; ++g)
          func(++idx, transport_view.MapDOF(i, m, first_group + g));
  }
}

} // namespace

template <typename Functor>
int64_t
LBSVecOps::GroupsetScopedCopy(LBSProblem& lbs_problem,
                              unsigned int gsi,
                              unsigned int gss,
                              Functor func)
{

  const auto& grid = lbs_problem.GetGrid();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto num_moments = lbs_problem.GetNumMoments();

  int64_t idx = -1;
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& transport_view = cell_transport_views[cell->local_id];
    const auto num_nodes = static_cast<std::size_t>(transport_view.GetNumNodes());
    for (std::size_t i = 0; i < num_nodes; ++i)
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

  auto& phi = (phi_opt == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()
                                                 : lbs_problem.GetPhiOldLocal();
  const auto& grid = lbs_problem.GetGrid();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  const auto num_groups = lbs_problem.GetNumGroups();

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& transport_view = cell_transport_views[cell->local_id];
    const auto num_nodes = static_cast<std::size_t>(transport_view.GetNumNodes());

    for (std::size_t i = 0; i < num_nodes; ++i)
      for (unsigned int g = 0; g < num_groups; ++g)
        phi[transport_view.MapDOF(i, 0, g)] = value;
  }
}

void
LBSVecOps::ScalePhiVector(LBSProblem& lbs_problem, PhiSTLOption phi_opt, double value)
{

  const auto& groupsets = lbs_problem.GetGroupsets();

  if (phi_opt == PhiSTLOption::PHI_NEW)
    lbs_problem.ScalePhiNew(value);
  else
    lbs_problem.ScalePhiOld(value);

  for (const auto& groupset : groupsets)
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
LBSVecOps::SetGSPETScVecFromPrimarySTLvector(
  LBSProblem& lbs_problem, const LBSGroupset& groupset, Vec dest, PhiSTLOption src, int64_t offset)
{
  const auto& src_phi =
    (src == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  PetscScalar* petsc_dest = nullptr;
  OpenSnPETScCall(VecGetArray(dest, &petsc_dest));
  int64_t index = offset + GroupsetScopedCopy(lbs_problem,
                                              groupset.first_group,
                                              groupset.GetNumGroups(),
                                              [&](int64_t idx, size_t mapped_idx)
                                              { petsc_dest[offset + idx] = src_phi[mapped_idx]; });
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
LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(
  LBSProblem& lbs_problem, const LBSGroupset& groupset, Vec src, PhiSTLOption dest, int64_t offset)
{
  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  const PetscScalar* petsc_src = nullptr;
  OpenSnPETScCall(VecGetArrayRead(src, &petsc_src));
  int64_t index = offset + GroupsetScopedCopy(lbs_problem,
                                              groupset.first_group,
                                              groupset.GetNumGroups(),
                                              [&](int64_t idx, size_t mapped_idx)
                                              { dest_phi[mapped_idx] = petsc_src[offset + idx]; });
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
  PetscScalar* petsc_dest = nullptr;
  OpenSnPETScCall(VecGetArray(dest, &petsc_dest));
  GroupRangeCopy(lbs_problem,
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
  const PetscScalar* petsc_src = nullptr;
  OpenSnPETScCall(VecGetArrayRead(src, &petsc_src));
  GroupRangeCopy(lbs_problem,
                 first_group_id,
                 last_group_id - first_group_id + 1,
                 [&](int64_t idx, size_t mapped_idx) { dest_phi[mapped_idx] = petsc_src[idx]; });
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
      groupset.angle_agg->SetDelayedPsiNew2Old();
    else if (src == PhiSTLOption::PHI_OLD and dest == PhiSTLOption::PHI_NEW)
      groupset.angle_agg->SetDelayedPsiOld2New();
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
  const auto& groupsets = lbs_problem.GetGroupsets();
  PetscScalar* x_ref = nullptr;
  OpenSnPETScCall(VecGetArray(x, &x_ref));

  int64_t offset = 0;
  for (auto gs_id : groupset_ids)
  {
    const auto& groupset = groupsets.at(gs_id);
    auto gsi = groupset.first_group;
    auto gsf = groupset.last_group;
    auto gss = gsf - gsi + 1;

    int64_t index = offset + GroupsetScopedCopy(lbs_problem,
                                                gsi,
                                                gss,
                                                [&](int64_t idx, size_t mapped_idx)
                                                { x_ref[offset + idx] = y[mapped_idx]; });
    if (groupset.angle_agg)
    {
      if (which_phi == PhiSTLOption::PHI_NEW)
        groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, x_ref);
      else if (which_phi == PhiSTLOption::PHI_OLD)
        groupset.angle_agg->AppendOldDelayedAngularDOFsToArray(index, x_ref);
    }
    offset = index + 1;
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
  const auto& groupsets = lbs_problem.GetGroupsets();
  const PetscScalar* x_ref = nullptr;
  OpenSnPETScCall(VecGetArrayRead(x, &x_ref));

  int64_t offset = 0;
  for (auto gs_id : groupset_ids)
  {
    const auto& groupset = groupsets.at(gs_id);
    auto gsi = groupset.first_group;
    auto gsf = groupset.last_group;
    auto gss = gsf - gsi + 1;

    int64_t index = offset + GroupsetScopedCopy(lbs_problem,
                                                gsi,
                                                gss,
                                                [&](int64_t idx, size_t mapped_idx)
                                                { y[mapped_idx] = x_ref[offset + idx]; });
    if (groupset.angle_agg)
    {
      if (which_phi == PhiSTLOption::PHI_NEW)
        groupset.angle_agg->SetNewDelayedAngularDOFsFromArray(index, x_ref);
      else if (which_phi == PhiSTLOption::PHI_OLD)
        groupset.angle_agg->SetOldDelayedAngularDOFsFromArray(index, x_ref);
    }
    offset = index + 1;
  }

  OpenSnPETScCall(VecRestoreArrayRead(x, &x_ref));
}

} // namespace opensn
