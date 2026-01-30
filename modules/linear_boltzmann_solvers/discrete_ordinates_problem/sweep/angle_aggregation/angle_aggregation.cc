// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

AngleAggregation::AngleAggregation(
  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  std::shared_ptr<AngularQuadrature>& quadrature,
  std::shared_ptr<MeshContinuum>& grid)
  : num_ang_unknowns_avail_(false), grid_(grid), quadrature_(quadrature), boundaries_(boundaries)
{
  for (const auto& bndry_id_cond : boundaries)
    bndry_id_cond.second->Setup(grid, *quadrature);
}

void
AngleAggregation::ZeroOutgoingDelayedPsi()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::ZeroOutgoingDelayedPsi");

  for (auto& angset : angle_set_groups_)
    for (auto& delayed_data : angset->GetFLUDS().DelayedPrelocIOutgoingPsi())
      std::fill(delayed_data.begin(), delayed_data.end(), 0.0);

  for (auto& angset : angle_set_groups_)
    std::fill(angset->GetFLUDS().DelayedLocalPsi().begin(),
              angset->GetFLUDS().DelayedLocalPsi().end(),
              0.0);
}

void
AngleAggregation::ZeroIncomingDelayedPsi()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::ZeroIncomingDelayedPsi");

  // Opposing reflecting bndries
  for (const auto& [bid, bndry] : boundaries_)
    bndry->ZeroOpposingDelayedAngularFluxOld();

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    std::fill(angle_set->GetFLUDS().DelayedLocalPsiOld().begin(),
              angle_set->GetFLUDS().DelayedLocalPsiOld().end(),
              0.0);

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
      std::fill(loc_vector.begin(), loc_vector.end(), 0.0);
}

void
AngleAggregation::InitializeReflectingBCs()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::InitializeReflectingBCs");

  const double epsilon = 1.0e-8;

  bool reflecting_bcs_initialized = false;

  for (auto& [bid, bndry] : boundaries_)
    bndry->InitializeDelayedAngularFlux(grid_, *quadrature_);

  for (auto& [bid, bndry] : boundaries_)
  {
    bndry->FinalizeDelayedAngularFluxSetup(bid, boundaries_);
    if (bndry->IsReflecting())
      reflecting_bcs_initialized = true;
  }

  if (reflecting_bcs_initialized)
    log.Log0Verbose1() << "Reflecting boundary conditions initialized.";
}

void
AngleAggregation::SetupAngleSetDependencies()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetupAngleSetDependencies");

  for (auto& angle_set : angle_set_groups_)
  {
    std::set<AngleSet*> following_angle_sets;
    for (const auto& [bid, bndry] : boundaries_)
      bndry->GetFollowingAngleSets(following_angle_sets, *this, *angle_set);
    angle_set->UpdateSweepDependencies(following_angle_sets);
  }
}

std::pair<size_t, size_t>
AngleAggregation::GetNumDelayedAngularDOFs()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::GetNumDelayedAngularDOFs");

  // Check if this is already developed
  if (num_ang_unknowns_avail_)
    return number_angular_unknowns_;

  // If not developed
  size_t local_ang_unknowns = 0;

  for (const auto& [bid, bndry] : boundaries_)
    local_ang_unknowns += bndry->CountDelayedAngularDOFsNew();

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    local_ang_unknowns += angle_set->GetFLUDS().DelayedLocalPsi().size();

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
      local_ang_unknowns += loc_vector.size();

  size_t global_ang_unknowns = 0;
  mpi_comm.all_reduce(local_ang_unknowns, global_ang_unknowns, mpi::op::sum<size_t>());

  number_angular_unknowns_ = {local_ang_unknowns, global_ang_unknowns};

  num_ang_unknowns_avail_ = true;
  return number_angular_unknowns_;
}

void
AngleAggregation::AppendNewDelayedAngularDOFsToArray(int64_t& index, double* x_ref)
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::AppendNewDelayedAngularDOFsToArray");

  for (auto& [bid, bndry] : boundaries_)
    bndry->AppendNewDelayedAngularDOFsToArray(index, x_ref);

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto val : angle_set->GetFLUDS().DelayedLocalPsi())
    {
      index++;
      x_ref[index] = val;
    }

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
      for (auto val : loc_vector)
      {
        index++;
        x_ref[index] = val;
      }
}

void
AngleAggregation::AppendOldDelayedAngularDOFsToArray(int64_t& index, double* x_ref)
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::AppendOldDelayedAngularDOFsToArray");

  for (auto& [bid, bndry] : boundaries_)
    bndry->AppendOldDelayedAngularDOFsToArray(index, x_ref);

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto val : angle_set->GetFLUDS().DelayedLocalPsiOld())
    {
      index++;
      x_ref[index] = val;
    }

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
      for (auto val : loc_vector)
      {
        index++;
        x_ref[index] = val;
      }
}

void
AngleAggregation::SetOldDelayedAngularDOFsFromArray(int64_t& index, const double* x_ref)
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetOldDelayedAngularDOFsFromArray");

  for (auto& [bid, bndry] : boundaries_)
    bndry->SetOldDelayedAngularDOFsFromArray(index, x_ref);

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& val : angle_set->GetFLUDS().DelayedLocalPsiOld())
    {
      index++;
      val = x_ref[index];
    }

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
      for (auto& val : loc_vector)
      {
        index++;
        val = x_ref[index];
      }
}

void
AngleAggregation::SetNewDelayedAngularDOFsFromArray(int64_t& index, const double* x_ref)
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetNewDelayedAngularDOFsFromArray");

  for (auto& [bid, bndry] : boundaries_)
    bndry->SetNewDelayedAngularDOFsFromArray(index, x_ref);

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& val : angle_set->GetFLUDS().DelayedLocalPsi())
    {
      index++;
      val = x_ref[index];
    }

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
      for (auto& val : loc_vector)
      {
        index++;
        val = x_ref[index];
      }
}

std::vector<double>
AngleAggregation::GetNewDelayedAngularDOFsAsSTLVector()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::GetNewDelayedAngularDOFsAsSTLVector");

  std::vector<double> psi_vector;

  auto psi_size = GetNumDelayedAngularDOFs();
  psi_vector.reserve(psi_size.first);

  for (auto& [bid, bndry] : boundaries_)
    bndry->AppendNewDelayedAngularDOFsToVector(psi_vector);

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto val : angle_set->GetFLUDS().DelayedLocalPsi())
      psi_vector.push_back(val);

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
      for (auto val : loc_vector)
        psi_vector.push_back(val);

  return psi_vector;
}

void
AngleAggregation::SetNewDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector)
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetNewDelayedAngularDOFsFromSTLVector");

  auto psi_size = GetNumDelayedAngularDOFs();
  size_t stl_size = stl_vector.size();
  if (stl_size != psi_size.first)
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": STL-vector size "
                           "is incompatible with number angular unknowns stored "
                           "in the angle-aggregation object.");

  size_t index = 0;
  for (auto& [bid, bndry] : boundaries_)
    bndry->SetNewDelayedAngularDOFsFromVector(stl_vector, index);

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& val : angle_set->GetFLUDS().DelayedLocalPsi())
      val = stl_vector[index++];

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
      for (auto& val : loc_vector)
        val = stl_vector[index++];
}

std::vector<double>
AngleAggregation::GetOldDelayedAngularDOFsAsSTLVector()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::GetOldDelayedAngularDOFsAsSTLVector");

  std::vector<double> psi_vector;

  auto psi_size = GetNumDelayedAngularDOFs();
  psi_vector.reserve(psi_size.first);

  for (auto& [bid, bndry] : boundaries_)
    bndry->AppendOldDelayedAngularDOFsToVector(psi_vector);

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto val : angle_set->GetFLUDS().DelayedLocalPsiOld())
      psi_vector.push_back(val);

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
      for (auto val : loc_vector)
        psi_vector.push_back(val);

  return psi_vector;
}

void
AngleAggregation::SetOldDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector)
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetOldDelayedAngularDOFsFromSTLVector");

  auto psi_size = GetNumDelayedAngularDOFs();
  size_t stl_size = stl_vector.size();
  if (stl_size != psi_size.first)
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": STL-vector size "
                           "is incompatible with number angular unknowns stored "
                           "in the angle-aggregation object.");

  size_t index = 0;
  for (auto& [bid, bndry] : boundaries_)
    bndry->SetOldDelayedAngularDOFsFromVector(stl_vector, index);

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& val : angle_set->GetFLUDS().DelayedLocalPsiOld())
      val = stl_vector[index++];

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
      for (auto& val : loc_vector)
        val = stl_vector[index++];
}

void
AngleAggregation::SetDelayedPsiOld2New()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetDelayedPsiOld2New");

  for (auto& [bid, bndry] : boundaries_)
    bndry->CopyDelayedAngularFluxOldToNew();

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    angle_set->GetFLUDS().SetDelayedLocalPsiOldToNew();

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    angle_set->GetFLUDS().SetDelayedOutgoingPsiOldToNew();
}

void
AngleAggregation::SetDelayedPsiNew2Old()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetDelayedPsiNew2Old");

  for (auto& [bid, bndry] : boundaries_)
    bndry->CopyDelayedAngularFluxNewToOld();

  // Intra-cell cycles
  for (auto& angle_set : angle_set_groups_)
    angle_set->GetFLUDS().SetDelayedLocalPsiNewToOld();

  // Inter location cycles
  for (auto& angle_set : angle_set_groups_)
    angle_set->GetFLUDS().SetDelayedOutgoingPsiNewToOld();
}

} // namespace opensn
