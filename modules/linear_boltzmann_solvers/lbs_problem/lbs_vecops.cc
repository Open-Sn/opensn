// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include <cstring>

namespace opensn
{

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
  CALI_CXX_MARK_SCOPE("LBSVecOps::SetGSPETScVecFromPrimarySTLvector");

  const auto& src_phi = (src == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()[groupset.id]
                                                       : lbs_problem.GetPhiOldLocal()[groupset.id];
  double* petsc_dest;
  VecGetArray(dest, &petsc_dest);
  std::memcpy(petsc_dest, src_phi.data(), src_phi.size() * sizeof(double));
  int64_t index = src_phi.size() - 1;
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
LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(LBSProblem& lbs_problem,
                                             const LBSGroupset& groupset,
                                             Vec src,
                                             PhiSTLOption dest)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::SetPrimarySTLvectorFromGSPETScVec");

  auto& dest_phi = (dest == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()[groupset.id]
                                                   : lbs_problem.GetPhiOldLocal()[groupset.id];
  const double* petsc_src;
  VecGetArrayRead(src, &petsc_src);
  std::memcpy(dest_phi.data(), petsc_src, dest_phi.size() * sizeof(double));
  int64_t index = dest_phi.size() - 1;
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
LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec(LBSProblem& lbs_problem,
                                                      const std::vector<LBSGroupset>& groupsets,
                                                      Vec src,
                                                      PhiSTLOption dest)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::SetPrimarySTLvectorFromGroupScopedPETScVec");

  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  const double* petsc_src;
  VecGetArrayRead(src, &petsc_src);
  int64_t ofst = 0;
  for (auto& gs : groupsets)
  {
    std::memcpy(dest_phi[gs.id].data(), petsc_src + ofst, dest_phi[gs.id].size() * sizeof(double));
    ofst += dest_phi[gs.id].size() * sizeof(double);
  }
  VecRestoreArrayRead(src, &petsc_src);
}

void
LBSVecOps::GSScopedCopyPrimarySTLvectors(LBSProblem& lbs_problem,
                                         const LBSGroupset& groupset,
                                         const std::vector<double>& src,
                                         std::vector<double>& dest)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::GSScopedCopyPrimarySTLvectors");
  assert(src.size() == dest.size());
  std::memcpy(dest.data(), src.data(), src.size() * sizeof(double));
}

void
LBSVecOps::GSScopedCopyPrimarySTLvectors(LBSProblem& lbs_problem,
                                         const LBSGroupset& groupset,
                                         PhiSTLOption src,
                                         PhiSTLOption dest)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::GSScopedCopyPrimarySTLvectors");

  const auto& src_phi =
    (src == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  auto& dest_phi =
    (dest == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal() : lbs_problem.GetPhiOldLocal();
  assert(dest_phi[groupset.id].size() == src_phi[groupset.id].size());
  std::memcpy(dest_phi[groupset.id].data(),
              src_phi[groupset.id].data(),
              src_phi[groupset.id].size() * sizeof(double));
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
                                                  const std::vector<int>& groupset_ids,
                                                  Vec x,
                                                  PhiSTLOption which_phi)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::SetMultiGSPETScVecFromPrimarySTLvector");

  auto& y = (which_phi == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()
                                                 : lbs_problem.GetPhiOldLocal();
  auto& groupsets = lbs_problem.GetGroupsets();
  double* x_ref;
  VecGetArray(x, &x_ref);

  int64_t index = 0;
  for (int gsi : groupset_ids)
  {
    std::memcpy(x_ref + index * sizeof(double), y[gsi].data(), y[gsi].size() * sizeof(double));
    index += y[gsi].size();
  }
  index--;

  for (int gsi : groupset_ids)
  {
    const auto& groupset = groupsets.at(gsi);
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
LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec(LBSProblem& lbs_problem,
                                                  const std::vector<int>& groupset_ids,
                                                  Vec x,
                                                  PhiSTLOption which_phi)
{
  CALI_CXX_MARK_SCOPE("LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec");

  auto& y = (which_phi == PhiSTLOption::PHI_NEW) ? lbs_problem.GetPhiNewLocal()
                                                 : lbs_problem.GetPhiOldLocal();
  auto& groupsets = lbs_problem.GetGroupsets();
  const double* x_ref;
  VecGetArrayRead(x, &x_ref);
  int64_t index = 0;
  for (int gsi : groupset_ids)
  {
    std::memcpy(y[gsi].data(), x_ref + index * sizeof(double), y[gsi].size() * sizeof(double));
    index += y[gsi].size();
  }
  index--;
  for (auto gsi : groupset_ids)
  {
    const auto& groupset = groupsets.at(gsi);
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
