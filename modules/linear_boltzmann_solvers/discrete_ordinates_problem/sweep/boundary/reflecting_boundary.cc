// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/reflecting_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

namespace opensn
{

template <bool ApplyOnOldFlux, typename Fn>
void
ReflectingBoundary::TransformGroupset(int groupset_id, Fn&& fn)
{
  if (not opposing_reflected_)
    return;
  auto&& f = std::forward<Fn>(fn);
  const auto& extra_data = extra_data_[groupset_id];
  const auto old_stride = extra_data.old_stride;
  double* flux = GetBoundaryFlux(groupset_id, ApplyOnOldFlux ? old_stride : 0);
  std::size_t size = old_stride * bank_[groupset_id].groupset_size;
  std::span<double> flux_view(flux, size);
  f(flux_view);
}

template <bool ApplyOnOldFlux, typename Fn>
void
ReflectingBoundary::TransformGroupset(int groupset_id, Fn&& fn) const
{
  if (not opposing_reflected_)
    return;
  auto&& f = std::forward<Fn>(fn);
  const auto& extra_data = extra_data_[groupset_id];
  const auto old_stride = extra_data.old_stride;
  const double* flux = GetBoundaryFlux(groupset_id, ApplyOnOldFlux ? old_stride : 0);
  std::size_t size = old_stride * bank_[groupset_id].groupset_size;
  std::span<const double> flux_view(flux, size);
  f(flux_view);
}

template <bool ApplyOnOldFlux, typename Fn>
void
ReflectingBoundary::ForEachDelayedAngularFlux(int groupset_id, Fn&& fn)
{
  if (not opposing_reflected_)
    return;
  auto&& f = std::forward<Fn>(fn);
  const auto& extra_data = extra_data_[groupset_id];
  const auto old_stride = extra_data.old_stride;
  double* flux = GetBoundaryFlux(groupset_id, ApplyOnOldFlux ? old_stride : 0);
  std::size_t size = old_stride * bank_[groupset_id].groupset_size;
  std::for_each_n(flux, size, f);
}

template <bool ApplyOnOldFlux, typename Fn>
void
ReflectingBoundary::ForEachDelayedAngularFlux(int groupset_id, Fn&& fn) const
{
  if (not opposing_reflected_)
    return;
  auto&& f = std::forward<Fn>(fn);
  const auto& extra_data = extra_data_[groupset_id];
  const auto old_stride = extra_data.old_stride;
  const double* flux = GetBoundaryFlux(groupset_id, ApplyOnOldFlux ? old_stride : 0);
  std::size_t size = old_stride * bank_[groupset_id].groupset_size;
  std::for_each_n(flux, size, f);
}

template <bool ApplyOnOldFlux, typename Fn>
void
ReflectingBoundary::ForEachDelayedAngularSlope(int groupset_id, Fn&& fn)
{
  if (not opposing_reflected_ or not delayed_angular_slope_enabled_)
    return;
  auto&& f = std::forward<Fn>(fn);
  const auto& extra_data = extra_data_[groupset_id];
  const auto old_stride = extra_data.old_stride;
  const auto offset = extra_data.slope_stride + (ApplyOnOldFlux ? old_stride : 0);
  double* slope = GetBoundaryFlux(groupset_id, offset);
  std::size_t size = old_stride * bank_[groupset_id].groupset_size;
  std::for_each_n(slope, size, f);
}

template <bool ApplyOnOldFlux, typename Fn>
void
ReflectingBoundary::ForEachDelayedAngularSlope(int groupset_id, Fn&& fn) const
{
  if (not opposing_reflected_ or not delayed_angular_slope_enabled_)
    return;
  auto&& f = std::forward<Fn>(fn);
  const auto& extra_data = extra_data_[groupset_id];
  const auto old_stride = extra_data.old_stride;
  const auto offset = extra_data.slope_stride + (ApplyOnOldFlux ? old_stride : 0);
  const double* slope = GetBoundaryFlux(groupset_id, offset);
  std::size_t size = old_stride * bank_[groupset_id].groupset_size;
  std::for_each_n(slope, size, f);
}

template <typename Fn>
void
ReflectingBoundary::AppendDelayedAngularFields(int groupset_id, bool old_store, Fn fn) const
{
  if (old_store)
  {
    ForEachDelayedAngularFlux<true>(groupset_id, fn);
    ForEachDelayedAngularSlope<true>(groupset_id, fn);
  }
  else
  {
    ForEachDelayedAngularFlux<false>(groupset_id, fn);
    ForEachDelayedAngularSlope<false>(groupset_id, fn);
  }
}

ReflectingBoundary::ReflectingBoundary(BoundaryBank& bank,
                                       std::uint64_t bid,
                                       const std::shared_ptr<MeshContinuum>& grid,
                                       const std::vector<LBSGroupset>& groupsets,
                                       const Vector3& normal,
                                       CoordinateSystemType coord_type)
  : SweepBoundary(bank, LBSBoundaryType::REFLECTING), normal_(normal), coord_type_(coord_type)
{
  extra_data_.reserve(groupsets.size());
  extra_data_.resize(groupsets.size());

  std::uint64_t face_node_counter = 0;
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_id = cell->local_id;
    for (unsigned int f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];
      if (not face.has_neighbor and face.neighbor_id == bid)
      {
        const auto num_face_nodes = face.vertex_ids.size();
        for (unsigned int fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          FaceNode fn(cell_id, f, fnode);
          facenode_to_index_[fn] = face_node_counter++;
        }
      }
    }
  }
}

void
ReflectingBoundary::InitializeReflectingMap(const std::vector<LBSGroupset>& groupsets)
{
  for (const auto& groupset : groupsets)
  {
    const auto& quadrature = groupset.quadrature;
    auto num_angles = quadrature->GetNumAngles();

    auto& extra_data = extra_data_[groupset.id];
    auto& reflected_anglenum = extra_data.reflected_anglenum;
    reflected_anglenum.reserve(num_angles);
    reflected_anglenum.resize(num_angles);

    const Vector3 ihat(1.0, 0.0, 0.0);
    const Vector3 jhat(0.0, 1.0, 0.0);

    for (std::uint32_t n = 0; n < num_angles; ++n)
    {
      const Vector3& omega_n = quadrature->GetOmega(n);
      Vector3 omega_reflected;

      switch (coord_type_)
      {
        case CoordinateSystemType::SPHERICAL:
          omega_reflected = -1.0 * omega_n;
          break;
        case CoordinateSystemType::CYLINDRICAL:
        {
          if (std::fabs(normal_.Dot(jhat)) > 0.999999 or normal_.Dot(ihat) < -0.999999)
            omega_reflected = omega_n - 2.0 * normal_ * omega_n.Dot(normal_);
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
      for (std::uint32_t nstar = 0; nstar < num_angles; ++nstar)
      {
        if (omega_reflected.Dot(quadrature->GetOmega(nstar)) > (1.0 - epsilon_))
        {
          reflected_anglenum[n] = nstar;
          found = true;
          break;
        }
      }

      if (not found)
      {
        throw std::logic_error("ReflectingBoundary: Reflected angle not found for angle " +
                               std::to_string(n) + " with direction " +
                               quadrature->GetOmega(n).PrintStr() +
                               ".\nThis can happen for two reasons:\n1) A quadrature is used "
                               "that is not symmetric about the axis associated with the "
                               "reflected boundary.\n2) The reflecting boundary is not "
                               "aligned with a reflecting axis of the quadrature.");
      }
    }
  }
}

void
ReflectingBoundary::InitializeAngleDependent(const std::vector<LBSGroupset>& groupsets)
{
  for (const auto& groupset : groupsets)
  {
    const auto& quadrature = groupset.quadrature;
    auto& angle_agg = *(groupset.angle_agg);
    auto num_angles = quadrature->GetNumAngles();
    auto& extra_data = extra_data_[groupset.id];
    auto& map_dirnum = extra_data.map_dirnum;
    auto& reflected_anglenum = extra_data.reflected_anglenum;
    auto& node_stride = extra_data.node_stride;
    auto& old_stride = extra_data.old_stride;
    auto& slope_stride = extra_data.slope_stride;

    map_dirnum.reserve(num_angles);
    map_dirnum.resize(num_angles, std::numeric_limits<std::uint64_t>::max());

    std::uint64_t internal_angle_idx = 0;
    for (const auto& angleset : angle_agg)
    {
      int inout_counter = 0;
      unsigned int orthogonal_counter = 0;
      for (const auto& angle : angleset->GetAngleIndices())
      {
        double dot = normal_.Dot(quadrature->GetOmega(angle));
        orthogonal_counter += (dot == 0);
        inout_counter += (dot > 0) - (dot < 0);
      }
      if (std::abs(inout_counter) + orthogonal_counter != angleset->GetAngleIndices().size())
        throw std::logic_error("ReflectingBoundary: Expect all angles in angle set to be all in or "
                               "all out.\n");

      bool is_outgoing = (inout_counter > 0);
      if (not is_outgoing)
        continue;
      for (const auto& angle : angleset->GetAngleIndices())
      {
        map_dirnum[angle] = internal_angle_idx;
        map_dirnum[reflected_anglenum[angle]] = internal_angle_idx;
        ++internal_angle_idx;
      }
    }

    auto& counter = bank_[groupset.id].counter;
    offset_[groupset.id] = counter;
    node_stride = internal_angle_idx;
    old_stride = facenode_to_index_.size() * node_stride;
    const std::uint64_t time_stride = opposing_reflected_ ? 2 * old_stride : old_stride;
    slope_stride = time_stride;
    const std::uint64_t total_stride =
      delayed_angular_slope_enabled_ ? 2 * time_stride : time_stride;
    counter += total_stride;
    bank_.ExtendBoundaryFlux(groupset.id, total_stride * groupset.GetNumGroups());
  }
}

void
ReflectingBoundary::GetFollowingAngleSets(int groupset_id,
                                          std::set<AngleSet*>& following_angle_sets,
                                          const AngleAggregation& angle_agg,
                                          const AngleSet& angleset)
{
  if (opposing_reflected_)
    return;

  const auto& omegas = angle_agg.GetQuadrature()->GetOmegas();
  for (const auto& angle_idx : angleset.GetAngleIndices())
  {
    if (omegas[angle_idx].Dot(normal_) < 0.0)
      continue;
    auto reflected_angle_num = extra_data_[groupset_id].reflected_anglenum[angle_idx];
    following_angle_sets.insert(angle_agg.GetAngleSetForAngleIndex(reflected_angle_num));
  }
}

void
ReflectingBoundary::SetOpposingReflected(
  std::uint64_t boundary_id,
  const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries)
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
}

void
ReflectingBoundary::ZeroOpposingDelayedAngularFluxOld(int groupset_id)
{
  TransformGroupset<true>(groupset_id, [](auto flux) { std::fill(flux.begin(), flux.end(), 0.0); });
  ForEachDelayedAngularSlope<true>(groupset_id, [](double& val) { val = 0.0; });
}

size_t
ReflectingBoundary::CountDelayedAngularDOFsNew(int groupset_id) const
{
  if (not opposing_reflected_)
    return 0;
  const auto& extra_data = extra_data_[groupset_id];
  const auto field_size = extra_data.old_stride * bank_[groupset_id].groupset_size;
  return delayed_angular_slope_enabled_ ? 2 * field_size : field_size;
}

size_t
ReflectingBoundary::CountDelayedAngularDOFsOld(int groupset_id) const
{
  return CountDelayedAngularDOFsNew(groupset_id);
}

void
ReflectingBoundary::AppendNewDelayedAngularDOFsToVector(int groupset_id,
                                                        std::vector<double>& output) const
{
  AppendDelayedAngularFields(groupset_id, false, [&output](double val) { output.push_back(val); });
}

void
ReflectingBoundary::AppendOldDelayedAngularDOFsToVector(int groupset_id,
                                                        std::vector<double>& output) const
{
  AppendDelayedAngularFields(groupset_id, true, [&output](double val) { output.push_back(val); });
}

void
ReflectingBoundary::AppendNewDelayedAngularDOFsToArray(int groupset_id,
                                                       int64_t& index,
                                                       double* buffer) const
{
  AppendDelayedAngularFields(
    groupset_id, false, [&index, buffer](double val) { buffer[++index] = val; });
}

void
ReflectingBoundary::AppendOldDelayedAngularDOFsToArray(int groupset_id,
                                                       int64_t& index,
                                                       double* buffer) const
{
  AppendDelayedAngularFields(
    groupset_id, true, [&index, buffer](double val) { buffer[++index] = val; });
}

void
ReflectingBoundary::SetNewDelayedAngularDOFsFromArray(int groupset_id,
                                                      int64_t& index,
                                                      const double* buffer)
{
  ForEachDelayedAngularFlux<false>(groupset_id,
                                   [&index, buffer](double& psi) { psi = buffer[++index]; });
  ForEachDelayedAngularSlope<false>(groupset_id,
                                    [&index, buffer](double& slope) { slope = buffer[++index]; });
}

void
ReflectingBoundary::SetOldDelayedAngularDOFsFromArray(int groupset_id,
                                                      int64_t& index,
                                                      const double* buffer)
{
  ForEachDelayedAngularFlux<true>(groupset_id,
                                  [&index, buffer](double& psi) { psi = buffer[++index]; });
  ForEachDelayedAngularSlope<true>(groupset_id,
                                   [&index, buffer](double& slope) { slope = buffer[++index]; });
}

void
ReflectingBoundary::SetNewDelayedAngularDOFsFromVector(int groupset_id,
                                                       const std::vector<double>& values,
                                                       size_t& index)
{
  ForEachDelayedAngularFlux<false>(groupset_id,
                                   [&values, &index](double& psi) { psi = values[index++]; });
  ForEachDelayedAngularSlope<false>(groupset_id,
                                    [&values, &index](double& slope) { slope = values[index++]; });
}

void
ReflectingBoundary::SetOldDelayedAngularDOFsFromVector(int groupset_id,
                                                       const std::vector<double>& values,
                                                       size_t& index)
{
  ForEachDelayedAngularFlux<true>(groupset_id,
                                  [&values, &index](double& psi) { psi = values[index++]; });
  ForEachDelayedAngularSlope<true>(groupset_id,
                                   [&values, &index](double& slope) { slope = values[index++]; });
}

void
ReflectingBoundary::CopyDelayedAngularFluxOldToNew(int groupset_id)
{
  if (not opposing_reflected_)
    return;
  const auto& extra_data = extra_data_[groupset_id];
  const auto field_size = extra_data.old_stride * bank_[groupset_id].groupset_size;
  double* flux = GetBoundaryFlux(groupset_id);
  const double* flux_old = flux + field_size;
  std::copy_n(flux_old, field_size, flux);
  if (delayed_angular_slope_enabled_)
  {
    double* slope = GetBoundaryFlux(groupset_id, extra_data.slope_stride);
    const double* slope_old = slope + field_size;
    std::copy_n(slope_old, field_size, slope);
  }
}

void
ReflectingBoundary::CopyDelayedAngularFluxNewToOld(int groupset_id)
{
  if (not opposing_reflected_)
    return;
  const auto& extra_data = extra_data_[groupset_id];
  const auto field_size = extra_data.old_stride * bank_[groupset_id].groupset_size;
  double* flux = GetBoundaryFlux(groupset_id);
  double* flux_old = flux + field_size;
  std::copy_n(flux, field_size, flux_old);
  if (delayed_angular_slope_enabled_)
  {
    double* slope = GetBoundaryFlux(groupset_id, extra_data.slope_stride);
    double* slope_old = slope + field_size;
    std::copy_n(slope, field_size, slope_old);
  }
}

double*
ReflectingBoundary::PsiIncoming(std::uint32_t cell_local_id,
                                unsigned int face_num,
                                unsigned int fi,
                                unsigned int angle_num,
                                int groupset_id,
                                unsigned int group_idx)
{
  FaceNode fn(cell_local_id, face_num, fi);
  std::uint64_t facenode_idx = facenode_to_index_.at(fn);
  auto& extra_data = extra_data_[groupset_id];
  const auto& node_stride = extra_data.node_stride;
  const auto& old_stride = extra_data.old_stride;
  std::uint64_t node_angle_offset = facenode_idx * node_stride + extra_data.map_dirnum[angle_num];
  if (opposing_reflected_)
    node_angle_offset += old_stride;
  return GetBoundaryFlux(groupset_id, node_angle_offset) + group_idx;
}

double*
ReflectingBoundary::PsiIncomingE(std::uint32_t cell_local_id,
                                 unsigned int face_num,
                                 unsigned int fi,
                                 unsigned int angle_num,
                                 int groupset_id,
                                 unsigned int group_num)
{
  if (not delayed_angular_slope_enabled_)
    return ZeroFlux(groupset_id, group_num);
  FaceNode fn(cell_local_id, face_num, fi);
  std::uint64_t facenode_idx = facenode_to_index_.at(fn);
  auto& extra_data = extra_data_[groupset_id];
  const auto& node_stride = extra_data.node_stride;
  const auto& old_stride = extra_data.old_stride;
  std::uint64_t node_angle_offset =
    extra_data.slope_stride + facenode_idx * node_stride + extra_data.map_dirnum[angle_num];
  if (opposing_reflected_)
    node_angle_offset += old_stride;
  return GetBoundaryFlux(groupset_id, node_angle_offset) + group_num;
}

double*
ReflectingBoundary::PsiOutgoing(uint64_t cell_local_id,
                                unsigned int face_num,
                                unsigned int fi,
                                unsigned int angle_num,
                                int groupset_id)
{
  FaceNode fn(cell_local_id, face_num, fi);
  std::uint64_t facenode_idx = facenode_to_index_.at(fn);
  auto& extra_data = extra_data_[groupset_id];
  const auto& node_stride = extra_data.node_stride;
  std::uint64_t node_angle_offset = facenode_idx * node_stride + extra_data.map_dirnum[angle_num];
  return GetBoundaryFlux(groupset_id, node_angle_offset);
}

double*
ReflectingBoundary::PsiOutgoingE(uint64_t cell_local_id,
                                 unsigned int face_num,
                                 unsigned int fi,
                                 unsigned int angle_num,
                                 int groupset_id)
{
  if (not delayed_angular_slope_enabled_)
    return ZeroFlux(groupset_id, 0);
  FaceNode fn(cell_local_id, face_num, fi);
  std::uint64_t facenode_idx = facenode_to_index_.at(fn);
  auto& extra_data = extra_data_[groupset_id];
  const auto& node_stride = extra_data.node_stride;
  std::uint64_t node_angle_offset =
    extra_data.slope_stride + facenode_idx * node_stride + extra_data.map_dirnum[angle_num];
  return GetBoundaryFlux(groupset_id, node_angle_offset);
}

std::uint64_t
ReflectingBoundary::GetOffsetToAngleset(const FaceNode& face_node,
                                        AngleSet& anglset,
                                        bool is_outgoing)
{
  std::uint64_t facenode_idx = facenode_to_index_.at(face_node);
  int groupset_id = anglset.GetGroupsetID();
  auto& extra_data = extra_data_[groupset_id];
  const auto& node_stride = extra_data.node_stride;
  const auto& old_stride = extra_data.old_stride;
  std::uint64_t node_angle_offset =
    facenode_idx * node_stride + extra_data.map_dirnum[anglset.GetAngleIndices()[0]];
  if (not is_outgoing and opposing_reflected_)
    node_angle_offset += old_stride;
  return (offset_[groupset_id] + node_angle_offset) * bank_[groupset_id].groupset_size;
}

} // namespace opensn
