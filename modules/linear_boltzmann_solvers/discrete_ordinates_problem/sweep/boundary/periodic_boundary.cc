// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/periodic_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace opensn
{

PeriodicBoundary::PeriodicBoundary(BoundaryBank& bank,
                                   std::uint64_t boundary_id,
                                   const std::shared_ptr<MeshContinuum>& grid,
                                   const std::vector<LBSGroupset>& groupsets)
  : SweepBoundary(bank, LBSBoundaryType::PERIODIC), boundary_id_(boundary_id)
{
  groupset_data_.resize(groupsets.size());
  for (const auto& cell : grid->GetLocalCells())
    for (unsigned int f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];
      if (face.has_neighbor or face.neighbor_id != boundary_id)
        continue;
      for (unsigned int fi = 0; fi < face.vertex_ids.size(); ++fi)
      {
        node_to_index_.emplace(FaceNode(cell->local_id, f, fi), node_to_index_.size());
        const auto& vertex = grid->GlobalVertex(face.vertex_ids[fi]);
        local_geometry_.insert(
          local_geometry_.end(),
          {vertex.x, vertex.y, vertex.z, face.centroid.x, face.centroid.y, face.centroid.z});
      }
    }
}

void
PeriodicBoundary::GatherNodeGeometry()
{
  const auto local_count = static_cast<std::uint64_t>(node_to_index_.size());
  std::vector<std::uint64_t> counts;
  mpi_comm.all_gather(local_count, counts);
  global_node_offsets_.assign(counts.size() + 1, 0);
  for (size_t rank = 0; rank < counts.size(); ++rank)
    global_node_offsets_[rank + 1] = global_node_offsets_[rank] + counts[rank];
  mpi_comm.all_gather(local_geometry_, global_geometry_);
  if (global_geometry_.size() != 6 * global_node_offsets_.back())
    throw std::logic_error("Periodic boundary geometry gather has an invalid size.");
}

void
PeriodicBoundary::MatchPartnerNodes()
{
  if (partner_ == nullptr or global_node_offsets_.empty() or partner_->global_node_offsets_.empty())
    throw std::logic_error("Periodic boundary pair was not initialized.");

  const auto node_count = global_node_offsets_.back();
  const auto partner_count = partner_->global_node_offsets_.back();
  if (node_count == 0 or node_count != partner_count)
    throw std::runtime_error("Periodic boundary " + std::to_string(boundary_id_) +
                             " has no one-to-one partner face nodes.");

  double translation[3] = {0.0, 0.0, 0.0};
  for (std::uint64_t i = 0; i < node_count; ++i)
    for (int d = 0; d < 3; ++d)
      translation[d] += partner_->global_geometry_[6 * i + 3 + d] - global_geometry_[6 * i + 3 + d];
  for (auto& value : translation)
    value /= static_cast<double>(node_count);

  double scale = 1.0;
  for (int d = 0; d < 3; ++d)
  {
    double lower = std::numeric_limits<double>::max();
    double upper = std::numeric_limits<double>::lowest();
    for (std::uint64_t i = 0; i < node_count; ++i)
      for (const auto* geometry : {&global_geometry_, &partner_->global_geometry_})
      {
        lower = std::min(lower, (*geometry)[6 * i + d]);
        upper = std::max(upper, (*geometry)[6 * i + d]);
      }
    scale = std::max(scale, upper - lower);
  }
  const double tolerance = 1.0e-9 * scale;

  std::vector<std::uint64_t> global_mapping(node_count, std::numeric_limits<std::uint64_t>::max());
  std::vector<std::uint64_t> own_order(node_count), partner_order(partner_count);
  std::iota(own_order.begin(), own_order.end(), 0);
  std::iota(partner_order.begin(), partner_order.end(), 0);
  const auto normal_axis = static_cast<int>(boundary_id_ / 2);
  const auto less_by_tangential_position =
    [normal_axis](const std::vector<double>& geometry, std::uint64_t a, std::uint64_t b)
  {
    // A face centroid distinguishes adjacent cells that share a boundary vertex.
    for (int offset : {3, 0})
      for (int d = 0; d < 3; ++d)
      {
        if (d == normal_axis)
          continue;
        const auto av = geometry[6 * a + offset + d];
        const auto bv = geometry[6 * b + offset + d];
        if (av != bv)
          return av < bv;
      }
    return false;
  };
  std::sort(own_order.begin(),
            own_order.end(),
            [&](auto a, auto b) { return less_by_tangential_position(global_geometry_, a, b); });
  std::sort(partner_order.begin(),
            partner_order.end(),
            [&](auto a, auto b)
            { return less_by_tangential_position(partner_->global_geometry_, a, b); });

  for (std::uint64_t position = 0; position < node_count; ++position)
  {
    const auto i = own_order[position];
    const auto j = partner_order[position];
    for (int d = 0; d < 3; ++d)
    {
      if (std::fabs(global_geometry_[6 * i + d] + translation[d] -
                    partner_->global_geometry_[6 * j + d]) > tolerance or
          std::fabs(global_geometry_[6 * i + 3 + d] + translation[d] -
                    partner_->global_geometry_[6 * j + 3 + d]) > tolerance)
        throw std::runtime_error("Periodic boundary faces do not match under translation.");
    }
    global_mapping[i] = j;
  }

  const auto rank = static_cast<size_t>(mpi_comm.rank());
  partner_node_index_.assign(global_mapping.begin() + global_node_offsets_[rank],
                             global_mapping.begin() + global_node_offsets_[rank + 1]);
}

void
PeriodicBoundary::InitializeAngleDependent(const std::vector<LBSGroupset>& groupsets)
{
  for (const auto& groupset : groupsets)
  {
    auto& data = groupset_data_[groupset.id];
    data.num_angles = groupset.quadrature->GetNumAngles();
    data.incoming_angle_index.assign(data.num_angles, std::numeric_limits<std::uint64_t>::max());
    const auto normal_axis = static_cast<int>(boundary_id_ / 2);
    const auto normal_sign = boundary_id_ % 2 == 0 ? -1.0 : 1.0;
    for (std::uint64_t angle = 0; angle < data.num_angles; ++angle)
    {
      const auto& omega = groupset.quadrature->GetOmega(angle);
      const double component = normal_axis == 0 ? omega.x : normal_axis == 1 ? omega.y : omega.z;
      if (normal_sign * component < 0.0)
        data.incoming_angle_index[angle] = data.num_incoming_angles++;
    }

    data.incoming_stride = node_to_index_.size() * data.num_incoming_angles;
    data.outgoing_stride = node_to_index_.size() * data.num_angles;
    auto& common = bank_[groupset.id];
    offset_[groupset.id] = common.counter;
    const auto slots = 2 * data.incoming_stride + data.outgoing_stride;
    common.counter += slots;
    bank_.ExtendBoundaryFlux(groupset.id, slots * common.groupset_size);
  }
}

size_t
PeriodicBoundary::IncomingSize(int groupset_id) const
{
  return groupset_data_[groupset_id].incoming_stride * bank_[groupset_id].groupset_size;
}

double*
PeriodicBoundary::IncomingNew(int groupset_id)
{
  return GetBoundaryFlux(groupset_id);
}

double*
PeriodicBoundary::IncomingOld(int groupset_id)
{
  return GetBoundaryFlux(groupset_id, groupset_data_[groupset_id].incoming_stride);
}

const double*
PeriodicBoundary::IncomingNew(int groupset_id) const
{
  return GetBoundaryFlux(groupset_id);
}

const double*
PeriodicBoundary::IncomingOld(int groupset_id) const
{
  return GetBoundaryFlux(groupset_id, groupset_data_[groupset_id].incoming_stride);
}

std::vector<double>
PeriodicBoundary::GatherOutgoing(int groupset_id) const
{
  const auto& data = groupset_data_[groupset_id];
  const auto group_count = bank_[groupset_id].groupset_size;
  const auto size = data.outgoing_stride * group_count;
  const auto* outgoing = GetBoundaryFlux(groupset_id, 2 * data.incoming_stride);
  std::vector<double> local(outgoing, outgoing + size);
  std::vector<double> global;
  mpi_comm.all_gather(local, global);
  if (global.size() != global_node_offsets_.back() * data.num_angles * group_count)
    throw std::logic_error("Periodic boundary outgoing flux gather has an invalid size.");
  return global;
}

void
PeriodicBoundary::AssembleIncoming(int groupset_id)
{
  const auto& data = groupset_data_[groupset_id];
  const auto group_count = bank_[groupset_id].groupset_size;
  const auto partner_outgoing = partner_->GatherOutgoing(groupset_id);
  auto* incoming = IncomingNew(groupset_id);
  for (size_t local_node = 0; local_node < partner_node_index_.size(); ++local_node)
    for (std::uint64_t angle = 0; angle < data.num_angles; ++angle)
    {
      const auto incoming_angle = data.incoming_angle_index[angle];
      if (incoming_angle == std::numeric_limits<std::uint64_t>::max())
        continue;
      const auto source = (partner_node_index_[local_node] * data.num_angles + angle) * group_count;
      const auto target = (local_node * data.num_incoming_angles + incoming_angle) * group_count;
      std::copy_n(partner_outgoing.data() + source, group_count, incoming + target);
    }
}

void
PeriodicBoundary::ZeroOpposingDelayedAngularFluxOld(int groupset_id)
{
  std::fill_n(IncomingOld(groupset_id), IncomingSize(groupset_id), 0.0);
}

size_t
PeriodicBoundary::CountDelayedAngularDOFsNew(int groupset_id) const
{
  return IncomingSize(groupset_id);
}

size_t
PeriodicBoundary::CountDelayedAngularDOFsOld(int groupset_id) const
{
  return IncomingSize(groupset_id);
}

void
PeriodicBoundary::PrepareNewDelayedAngularFlux(int groupset_id)
{
  AssembleIncoming(groupset_id);
}

void
PeriodicBoundary::AppendNewDelayedAngularDOFsToVector(int groupset_id,
                                                      std::vector<double>& output) const
{
  const auto* values = IncomingNew(groupset_id);
  output.insert(output.end(), values, values + IncomingSize(groupset_id));
}

void
PeriodicBoundary::AppendOldDelayedAngularDOFsToVector(int groupset_id,
                                                      std::vector<double>& output) const
{
  const auto* values = IncomingOld(groupset_id);
  output.insert(output.end(), values, values + IncomingSize(groupset_id));
}

void
PeriodicBoundary::AppendNewDelayedAngularDOFsToArray(int groupset_id,
                                                     int64_t& index,
                                                     double* buffer) const
{
  for (size_t i = 0; i < IncomingSize(groupset_id); ++i)
    buffer[++index] = IncomingNew(groupset_id)[i];
}

void
PeriodicBoundary::AppendOldDelayedAngularDOFsToArray(int groupset_id,
                                                     int64_t& index,
                                                     double* buffer) const
{
  for (size_t i = 0; i < IncomingSize(groupset_id); ++i)
    buffer[++index] = IncomingOld(groupset_id)[i];
}

void
PeriodicBoundary::SetNewDelayedAngularDOFsFromArray(int groupset_id,
                                                    int64_t& index,
                                                    const double* buffer)
{
  for (size_t i = 0; i < IncomingSize(groupset_id); ++i)
    IncomingNew(groupset_id)[i] = buffer[++index];
}

void
PeriodicBoundary::SetOldDelayedAngularDOFsFromArray(int groupset_id,
                                                    int64_t& index,
                                                    const double* buffer)
{
  for (size_t i = 0; i < IncomingSize(groupset_id); ++i)
    IncomingOld(groupset_id)[i] = buffer[++index];
}

void
PeriodicBoundary::SetNewDelayedAngularDOFsFromVector(int groupset_id,
                                                     const std::vector<double>& values,
                                                     size_t& index)
{
  for (size_t i = 0; i < IncomingSize(groupset_id); ++i)
    IncomingNew(groupset_id)[i] = values[index++];
}

void
PeriodicBoundary::SetOldDelayedAngularDOFsFromVector(int groupset_id,
                                                     const std::vector<double>& values,
                                                     size_t& index)
{
  for (size_t i = 0; i < IncomingSize(groupset_id); ++i)
    IncomingOld(groupset_id)[i] = values[index++];
}

void
PeriodicBoundary::CopyDelayedAngularFluxOldToNew(int groupset_id)
{
  std::copy_n(IncomingOld(groupset_id), IncomingSize(groupset_id), IncomingNew(groupset_id));
}

void
PeriodicBoundary::CopyDelayedAngularFluxNewToOld(int groupset_id)
{
  std::copy_n(IncomingNew(groupset_id), IncomingSize(groupset_id), IncomingOld(groupset_id));
}

double*
PeriodicBoundary::PsiIncoming(std::uint32_t cell_local_id,
                              unsigned int face_num,
                              unsigned int fi,
                              unsigned int angle_num,
                              int groupset_id,
                              unsigned int group_idx)
{
  const auto& data = groupset_data_[groupset_id];
  const auto node = node_to_index_.at(FaceNode(cell_local_id, face_num, fi));
  const auto incoming_angle = data.incoming_angle_index.at(angle_num);
  if (incoming_angle == std::numeric_limits<std::uint64_t>::max())
    throw std::logic_error("Periodic boundary queried for a non-incoming angle.");
  return IncomingOld(groupset_id) +
         (node * data.num_incoming_angles + incoming_angle) * bank_[groupset_id].groupset_size +
         group_idx;
}

double*
PeriodicBoundary::PsiOutgoing(uint64_t cell_local_id,
                              unsigned int face_num,
                              unsigned int fi,
                              unsigned int angle_num,
                              int groupset_id)
{
  const auto& data = groupset_data_[groupset_id];
  const auto node = node_to_index_.at(FaceNode(cell_local_id, face_num, fi));
  return GetBoundaryFlux(groupset_id,
                         2 * data.incoming_stride + node * data.num_angles + angle_num);
}

std::uint64_t
PeriodicBoundary::GetOffsetToAngleset(const FaceNode& face_node,
                                      AngleSet& angleset,
                                      bool is_outgoing)
{
  const auto groupset_id = angleset.GetGroupsetID();
  const auto& data = groupset_data_[groupset_id];
  const auto node = node_to_index_.at(face_node);
  const auto angle = angleset.GetAngleIndices().front();
  const auto slot = is_outgoing ? 2 * data.incoming_stride + node * data.num_angles + angle
                                : data.incoming_stride + node * data.num_incoming_angles +
                                    data.incoming_angle_index.at(angle);
  return (offset_[groupset_id] + slot) * bank_[groupset_id].groupset_size;
}

} // namespace opensn
