// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/boundary/reflecting_boundary.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{
namespace lbs
{

AngleAggregation::AngleAggregation(
  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  size_t num_groups,
  size_t num_group_subsets,
  std::shared_ptr<AngularQuadrature>& quadrature,
  std::shared_ptr<MeshContinuum>& grid)
  : is_setup_(false),
    num_groups_(num_groups),
    num_group_subsets_(num_group_subsets),
    num_ang_unknowns_avail_(false),
    grid_(grid),
    quadrature_(quadrature),
    boundaries_(boundaries)
{
  for (auto& bndry_id_cond : boundaries)
    bndry_id_cond.second->Setup(*grid, *quadrature);

  is_setup_ = true;
}

void
AngleAggregation::ZeroOutgoingDelayedPsi()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::ZeroOutgoingDelayedPsi");

  for (auto& angsetgrp : angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      for (auto& delayed_data : angset->GetFLUDS().DelayedPrelocIOutgoingPsi())
        Set(delayed_data, 0.0);

  for (auto& angsetgrp : angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      Set(angset->GetFLUDS().DelayedLocalPsi(), 0.0);
}

void
AngleAggregation::ZeroIncomingDelayedPsi()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::ZeroIncomingDelayedPsi");

  // Opposing reflecting bndries
  for (const auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                  val = 0.0;

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      Set(angle_set->GetFLUDS().DelayedLocalPsiOld(), 0.0);

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
        Set(loc_vector, 0.0);
}

void
AngleAggregation::InitializeReflectingBCs()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::InitializeReflectingBCs");

  const std::string fname = "AngleAggregation";
  const double epsilon = 1.0e-8;

  bool reflecting_bcs_initialized = false;

  const Vector3 ihat(1.0, 0.0, 0.0);
  const Vector3 jhat(0.0, 1.0, 0.0);

  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      size_t tot_num_angles = quadrature_->abscissae_.size();
      size_t num_local_cells = grid_->local_cells.size();
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      const auto& normal = rbndry.Normal();

      rbndry.GetReflectedAngleIndexMap().resize(tot_num_angles, -1);
      rbndry.GetAngleReadyFlags().resize(tot_num_angles,
                                         std::vector<bool>(num_group_subsets_, false));

      // Determine reflected angle and check that it is within the quadrature
      typedef Vector3 Vec3;
      for (int n = 0; n < tot_num_angles; ++n)
      {
        const Vec3& omega_n = quadrature_->omegas_[n];
        Vec3 omega_reflected;

        switch (rbndry.CoordType())
        {
          case CoordinateSystemType::SPHERICAL:
            omega_reflected = -1.0 * omega_n;
            break;
          case CoordinateSystemType::CYLINDRICAL:
          {
            // left, top and bottom is regular reflecting
            if (std::fabs(normal.Dot(jhat)) > 0.999999 or normal.Dot(ihat) < -0.999999)
              omega_reflected = omega_n - 2.0 * normal * omega_n.Dot(normal);
            // right derive their normal from omega_n
            else if (normal.Dot(ihat) > 0.999999)
            {
              Vec3 normal_star;
              if (omega_n.Dot(normal) > 0.0)
                normal_star = Vec3(omega_n.x, 0.0, omega_n.z).Normalized();
              else
                normal_star = Vec3(-omega_n.x, 0.0, -omega_n.y).Normalized();

              omega_reflected = omega_n - 2.0 * normal_star * omega_n.Dot(normal_star);
            }
          }
          break;
          case CoordinateSystemType::CARTESIAN:
          default:
            omega_reflected = omega_n - 2.0 * normal * omega_n.Dot(normal);
            break;
        }

        auto& index_map = rbndry.GetReflectedAngleIndexMap();
        for (int nstar = 0; nstar < tot_num_angles; ++nstar)
          if (omega_reflected.Dot(quadrature_->omegas_[nstar]) > (1.0 - epsilon))
          {
            index_map[n] = nstar;
            break;
          }

        if (index_map[n] < 0)
          throw std::logic_error(
            fname + ": Reflected angle not found for angle " + std::to_string(n) +
            " with direction " + quadrature_->omegas_[n].PrintStr() +
            ". This can happen for two reasons: i) A quadrature is used that is not symmetric "
            "about the axis associated with the reflected boundary, or ii) the reflecting boundary "
            "is not aligned with any reflecting axis of the quadrature.");
      }

      // Initialize storage for all outbound directions
      auto& heteroflux_new = rbndry.GetBoundaryFluxNew();
      auto& heteroflux_old = rbndry.GetBoundaryFluxOld();
      heteroflux_new.clear();
      heteroflux_old.clear();
      heteroflux_new.resize(tot_num_angles);
      for (int n = 0; n < tot_num_angles; ++n)
      {
        // Only continue if omega is outgoing
        if (quadrature_->omegas_[n].Dot(rbndry.Normal()) < 0.0)
          continue;

        // For cells
        auto& cell_vec = heteroflux_new[n];
        cell_vec.resize(num_local_cells);
        for (const auto& cell : grid_->local_cells)
        {
          const uint64_t c = cell.local_id_;

          // Check cell on ref bndry
          bool on_ref_bndry = false;
          for (const auto& face : cell.faces_)
          {
            if ((not face.has_neighbor_) and (face.normal_.Dot(rbndry.Normal()) > 0.999999))
            {
              on_ref_bndry = true;
              break;
            }
          }
          if (not on_ref_bndry)
            continue;

          // If cell on ref bndry
          cell_vec[c].resize(cell.faces_.size());
          int f = 0;
          for (const auto& face : cell.faces_)
          {
            if ((not face.has_neighbor_) and (face.normal_.Dot(rbndry.Normal()) > 0.999999))
            {
              cell_vec[c][f].clear();
              cell_vec[c][f].resize(face.vertex_ids_.size(), std::vector<double>(num_groups_, 0.0));
            }
            ++f;
          }
        } // for cells
      }   // for angles

      // Determine if boundary is opposing reflecting
      // The boundary with the smallest bid will
      // be marked as "opposing-reflecting" while
      // the other one will be just a regular
      // reflecting boundary
      for (const auto& [otherbid, otherbndry] : boundaries_)
      {
        if (bid == otherbid)
          continue;
        if (not otherbndry->IsReflecting())
          continue;

        const auto& otherRbndry = dynamic_cast<const ReflectingBoundary&>(*otherbndry);

        if (rbndry.Normal().Dot(otherRbndry.Normal()) < (0.0 - epsilon))
          if (bid < otherbid)
            rbndry.SetOpposingReflected(true);
      }

      if (rbndry.IsOpposingReflected())
        rbndry.GetBoundaryFluxOld() = rbndry.GetBoundaryFluxNew();

      reflecting_bcs_initialized = true;
    } // if reflecting
  }   // for bndry

  if (reflecting_bcs_initialized)
    log.Log0Verbose1() << "Reflecting boundary conditions initialized.";
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

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                local_ang_unknowns += dofvec.size();

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      local_ang_unknowns += angle_set->GetFLUDS().DelayedLocalPsi().size();

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
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

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                {
                  index++;
                  x_ref[index] = val;
                }

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto val : angle_set->GetFLUDS().DelayedLocalPsi())
      {
        index++;
        x_ref[index] = val;
      }

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
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

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                {
                  index++;
                  x_ref[index] = val;
                }

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto val : angle_set->GetFLUDS().DelayedLocalPsiOld())
      {
        index++;
        x_ref[index] = val;
      }

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
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

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                {
                  index++;
                  val = x_ref[index];
                }

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& val : angle_set->GetFLUDS().DelayedLocalPsiOld())
      {
        index++;
        val = x_ref[index];
      }

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
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

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                {
                  index++;
                  val = x_ref[index];
                }

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& val : angle_set->GetFLUDS().DelayedLocalPsi())
      {
        index++;
        val = x_ref[index];
      }

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
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

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                  psi_vector.push_back(val);

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto val : angle_set->GetFLUDS().DelayedLocalPsi())
        psi_vector.push_back(val);

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
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
  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                  val = stl_vector[index++];

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& val : angle_set->GetFLUDS().DelayedLocalPsi())
        val = stl_vector[index++];

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
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

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                  psi_vector.push_back(val);

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto val : angle_set->GetFLUDS().DelayedLocalPsiOld())
        psi_vector.push_back(val);

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
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
  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                  val = stl_vector[index++];

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& val : angle_set->GetFLUDS().DelayedLocalPsiOld())
        val = stl_vector[index++];

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
        for (auto& val : loc_vector)
          val = stl_vector[index++];
}

void
AngleAggregation::SetDelayedPsiOld2New()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetDelayedPsiOld2New");

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        rbndry.GetBoundaryFluxNew() = rbndry.GetBoundaryFluxOld();

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      angle_set->GetFLUDS().DelayedLocalPsi() = angle_set->GetFLUDS().DelayedLocalPsiOld();

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi() =
        angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld();
}

void
AngleAggregation::SetDelayedPsiNew2Old()
{
  CALI_CXX_MARK_SCOPE("AngleAggregation::SetDelayedPsiNew2Old");

  // Opposing reflecting bndries
  for (auto& [bid, bndry] : boundaries_)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (ReflectingBoundary&)(*bndry);

      if (rbndry.IsOpposingReflected())
        rbndry.GetBoundaryFluxOld() = rbndry.GetBoundaryFluxNew();

    } // if reflecting
  }   // for bndry

  // Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      angle_set->GetFLUDS().DelayedLocalPsiOld() = angle_set->GetFLUDS().DelayedLocalPsi();

  // Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld() =
        angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi();
}

} // namespace lbs
} // namespace opensn
