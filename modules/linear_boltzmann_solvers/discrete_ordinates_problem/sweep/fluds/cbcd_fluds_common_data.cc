// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"

namespace opensn
{

CBCD_FLUDSCommonData::CBCD_FLUDSCommonData(
  const SPDS& spds,
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SpatialDiscretization& sdm)
  : FLUDSCommonData(spds, grid_nodal_mappings)
{
  CopyFlattenNodeIndexToDevice(sdm);
}

CBCD_FLUDSCommonData::~CBCD_FLUDSCommonData()
{
  DeallocateDeviceMemory();
}

#ifndef __OPENSN_WITH_GPU__
void
CBCD_FLUDSCommonData::CopyFlattenNodeIndexToDevice(const SpatialDiscretization& sdm)
{
}

void
CBCD_FLUDSCommonData::DeallocateDeviceMemory()
{
}
#endif

const CBCD_IncomingNonlocalFace&
CBCD_FLUDSCommonData::GetIncomingNonlocalFace(const std::uint32_t source_partition_index,
                                              const std::uint32_t incoming_face_index) const
{
  const auto begin = source_to_incoming_face_offsets_[source_partition_index];
  return incoming_nonlocal_faces_[incoming_face_indices_by_source_[begin + incoming_face_index]];
}

} // namespace opensn
