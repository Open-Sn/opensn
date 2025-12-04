// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/reflecting_boundary.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace opensn
{

template <typename Fn>
void
ReflectingBoundary::ForEachDelayedAngularFlux(bool use_old_store, Fn&& fn)
{
  auto&& apply = std::forward<Fn>(fn);
  auto& flux = use_old_store ? boundary_flux_old_ : boundary_flux_;
  for (auto& angle : flux)
    for (auto& cellvec : angle)
      for (auto& facevec : cellvec)
        for (auto& dofvec : facevec)
          for (auto& val : dofvec)
            apply(val);
}

template <typename Fn>
void
ReflectingBoundary::ForEachDelayedAngularFluxConst(bool use_old_store, Fn&& fn) const
{
  auto&& apply = std::forward<Fn>(fn);
  const auto& flux = use_old_store ? boundary_flux_old_ : boundary_flux_;
  for (const auto& angle : flux)
    for (const auto& cellvec : angle)
      for (const auto& facevec : cellvec)
        for (const auto& dofvec : facevec)
          for (const auto& val : dofvec)
            apply(val);
}

void
ReflectingBoundary::InitializeDelayedAngularFlux(const std::shared_ptr<MeshContinuum>& grid,
                                                 const AngularQuadrature& quadrature)
{
  if (not grid)
    throw std::logic_error("ReflectingBoundary: mesh continuum is null.");

  const auto tot_num_angles = static_cast<int>(quadrature.abscissae.size());
  reflected_anglenum_.assign(tot_num_angles, -1);
  angle_readyflags_.assign(tot_num_angles, false);

  const Vector3 ihat(1.0, 0.0, 0.0);
  const Vector3 jhat(0.0, 1.0, 0.0);

  for (int n = 0; n < tot_num_angles; ++n)
  {
    const Vector3& omega_n = quadrature.omegas[n];
    Vector3 omega_reflected;

    switch (GetCoordType())
    {
      case CoordinateSystemType::SPHERICAL:
        omega_reflected = -1.0 * omega_n;
        break;
      case CoordinateSystemType::CYLINDRICAL:
      {
        // left, top and bottom are regular reflecting
        if (std::fabs(normal_.Dot(jhat)) > 0.999999 or normal_.Dot(ihat) < -0.999999)
          omega_reflected = omega_n - 2.0 * normal_ * omega_n.Dot(normal_);
        // right derives normal from omega_n
        else if (normal_.Dot(ihat) > 0.999999)
        {
          Vector3 normal_star;
          if (omega_n.Dot(normal_) > 0.0)
            normal_star = Vector3(omega_n.x, 0.0, omega_n.z).Normalized();
          else
            normal_star = Vector3(-omega_n.x, 0.0, -omega_n.y).Normalized();
          omega_reflected = omega_n - 2.0 * normal_star * omega_n.Dot(normal_star);
        }
        break;
      }
      case CoordinateSystemType::CARTESIAN:
      default:
        omega_reflected = omega_n - 2.0 * normal_ * omega_n.Dot(normal_);
        break;
    }

    bool found = false;
    for (int nstar = 0; nstar < tot_num_angles; ++nstar)
    {
      if (omega_reflected.Dot(quadrature.omegas[nstar]) > (1.0 - epsilon_))
      {
        reflected_anglenum_[n] = nstar;
        found = true;
        break;
      }
    }

    if (not found)
    {
      throw std::logic_error("ReflectingBoundary: Reflected angle not found for angle " +
                             std::to_string(n) + " with direction " +
                             quadrature.omegas[n].PrintStr() +
                             ".\nThis can happen for two reasons:\n1) A quadrature is used "
                             "that is not symmetric about the axis associated with the "
                             "reflected boundary.\n2) The reflecting boundary is not "
                             "aligned with a reflecting axis of the quadrature.");
    }
  }

  boundary_flux_.clear();
  boundary_flux_old_.clear();
  boundary_flux_.resize(tot_num_angles);

  const size_t num_local_cells = grid->local_cells.size();
  for (int n = 0; n < tot_num_angles; ++n)
  {
    if (quadrature.omegas[n].Dot(normal_) < 0.0)
      continue;

    auto& cell_vec = boundary_flux_[n];
    cell_vec.resize(num_local_cells);

    for (const auto& cell : grid->local_cells)
    {
      const uint64_t c = cell.local_id;
      bool on_ref_bndry = false;

      for (const auto& face : cell.faces)
      {
        if ((not face.has_neighbor) and (face.normal.Dot(normal_) > 0.999999))
        {
          on_ref_bndry = true;
          break;
        }
      }

      if (not on_ref_bndry)
        continue;

      cell_vec[c].resize(cell.faces.size());
      int f = 0;
      for (const auto& face : cell.faces)
      {
        if ((not face.has_neighbor) and (face.normal.Dot(normal_) > 0.999999))
        {
          cell_vec[c][f].clear();
          cell_vec[c][f].resize(face.vertex_ids.size(), std::vector<double>(num_groups_, 0.0));
        }
        ++f;
      }
    }
  }
}

void
ReflectingBoundary::FinalizeDelayedAngularFluxSetup(
  uint64_t boundary_id, const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries)
{
  for (const auto& [otherbid, otherbndry] : boundaries)
  {
    if (boundary_id == otherbid)
      continue;

    const Vector3* other_normal = otherbndry->GetNormalForReflection();
    if (other_normal == nullptr)
      continue;

    if (normal_.Dot(*other_normal) < (0.0 - epsilon_))
    {
      if (boundary_id < otherbid)
      {
        opposing_reflected_ = true;
        break;
      }
    }
  }

  if (opposing_reflected_)
    boundary_flux_old_ = boundary_flux_;
}

void
ReflectingBoundary::ZeroOpposingDelayedAngularFluxOld()
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFlux(true, [](double& val) { val = 0.0; });
}

size_t
ReflectingBoundary::CountDelayedAngularDOFsNew() const
{
  if (not opposing_reflected_)
    return 0;

  size_t counter = 0;
  ForEachDelayedAngularFluxConst(false, [&](const double&) { ++counter; });
  return counter;
}

size_t
ReflectingBoundary::CountDelayedAngularDOFsOld() const
{
  if (not opposing_reflected_)
    return 0;

  size_t counter = 0;
  ForEachDelayedAngularFluxConst(true, [&](const double&) { ++counter; });
  return counter;
}

void
ReflectingBoundary::AppendNewDelayedAngularDOFsToVector(std::vector<double>& output) const
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFluxConst(false, [&](const double& val) { output.push_back(val); });
}

void
ReflectingBoundary::AppendOldDelayedAngularDOFsToVector(std::vector<double>& output) const
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFluxConst(true, [&](const double& val) { output.push_back(val); });
}

void
ReflectingBoundary::AppendNewDelayedAngularDOFsToArray(int64_t& index, double* buffer) const
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFluxConst(false,
                                 [&](const double& val)
                                 {
                                   ++index;
                                   buffer[index] = val;
                                 });
}

void
ReflectingBoundary::AppendOldDelayedAngularDOFsToArray(int64_t& index, double* buffer) const
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFluxConst(true,
                                 [&](const double& val)
                                 {
                                   ++index;
                                   buffer[index] = val;
                                 });
}

void
ReflectingBoundary::SetNewDelayedAngularDOFsFromArray(int64_t& index, const double* buffer)
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFlux(false,
                            [&](double& val)
                            {
                              ++index;
                              val = buffer[index];
                            });
}

void
ReflectingBoundary::SetOldDelayedAngularDOFsFromArray(int64_t& index, const double* buffer)
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFlux(true,
                            [&](double& val)
                            {
                              ++index;
                              val = buffer[index];
                            });
}

void
ReflectingBoundary::SetNewDelayedAngularDOFsFromVector(const std::vector<double>& values,
                                                       size_t& index)
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFlux(false, [&](double& val) { val = values[index++]; });
}

void
ReflectingBoundary::SetOldDelayedAngularDOFsFromVector(const std::vector<double>& values,
                                                       size_t& index)
{
  if (not opposing_reflected_)
    return;

  ForEachDelayedAngularFlux(true, [&](double& val) { val = values[index++]; });
}

void
ReflectingBoundary::CopyDelayedAngularFluxOldToNew()
{
  if (not opposing_reflected_)
    return;

  boundary_flux_ = boundary_flux_old_;
}

void
ReflectingBoundary::CopyDelayedAngularFluxNewToOld()
{
  if (not opposing_reflected_)
    return;

  boundary_flux_old_ = boundary_flux_;
}

double*
ReflectingBoundary::PsiIncoming(uint64_t cell_local_id,
                                unsigned int face_num,
                                unsigned int fi,
                                unsigned int angle_num,
                                int group_num)
{
  int reflected_angle_num = reflected_anglenum_[angle_num];

  if (opposing_reflected_)
    return &boundary_flux_old_[reflected_angle_num][cell_local_id][face_num][fi].front();

  return &boundary_flux_[reflected_angle_num][cell_local_id][face_num][fi].front();
}

double*
ReflectingBoundary::PsiOutgoing(uint64_t cell_local_id,
                                unsigned int face_num,
                                unsigned int fi,
                                unsigned int angle_num)
{
  return &boundary_flux_[angle_num][cell_local_id][face_num][fi].front();
}

void
ReflectingBoundary::UpdateAnglesReadyStatus(const std::vector<std::uint32_t>& angles)
{
  for (const size_t n : angles)
    angle_readyflags_[reflected_anglenum_[n]] = true;
}

bool
ReflectingBoundary::CheckAnglesReadyStatus(const std::vector<std::uint32_t>& angles)
{
  if (opposing_reflected_)
    return true;

  for (const auto& n : angles)
    if (not boundary_flux_[reflected_anglenum_[n]].empty())
      if (not angle_readyflags_[n])
        return false;

  return true;
}

void
ReflectingBoundary::ResetAnglesReadyStatus()
{
  boundary_flux_old_ = boundary_flux_;
  std::fill(angle_readyflags_.begin(), angle_readyflags_.end(), false);
}

} // namespace opensn
