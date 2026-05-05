// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "framework/utils/error.h"
#include <cassert>
#include <algorithm>

namespace opensn
{

class SpatialDiscretization;

CBCD_FLUDSCommonData::CBCD_FLUDSCommonData(
  const SPDS& spds,
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SpatialDiscretization& sdm)
  : FLUDSCommonData(spds, grid_nodal_mappings),
    num_incoming_boundary_nodes_(0),
    num_outgoing_boundary_nodes_(0),
    num_incoming_nonlocal_faces_(0),
    num_incoming_nonlocal_nodes_(0),
    num_outgoing_nonlocal_faces_(0),
    num_outgoing_nonlocal_nodes_(0),
    device_cell_face_node_map_(nullptr)
{
  CopyFlattenedNodeIndexToDevice(sdm);
}

CBCD_FLUDSCommonData::~CBCD_FLUDSCommonData()
{
  DeallocateDeviceMemory();
}

#ifndef __OPENSN_WITH_GPU__
void
CBCD_FLUDSCommonData::CopyFlattenedNodeIndexToDevice(const SpatialDiscretization& sdm)
{
}

void
CBCD_FLUDSCommonData::DeallocateDeviceMemory()
{
}
#endif

const GroupedIncomingNonlocalFace&
CBCD_FLUDSCommonData::GetIncomingNonlocalFace(const std::uint32_t source_slot,
                                              const std::uint32_t source_face_index) const
{
  const auto begin = source_to_incoming_face_offsets_[source_slot];
  assert(begin + source_face_index < source_to_incoming_face_offsets_[source_slot + 1]);
  return incoming_nonlocal_faces_[incoming_face_indices_by_source_[begin + source_face_index]];
}

} // namespace opensn
