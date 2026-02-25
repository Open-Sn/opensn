// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"

namespace opensn
{

class SpatialDiscretization;

CBCD_FLUDSCommonData::CBCD_FLUDSCommonData(
  const SPDS& spds,
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SpatialDiscretization& sdm)
  : CBC_FLUDSCommonData(spds, grid_nodal_mappings),
    num_incoming_boundary_nodes_(0),
    num_outgoing_boundary_nodes_(0),
    num_incoming_nonlocal_nodes_(0),
    num_outgoing_nonlocal_nodes_(0),
    device_cell_face_node_map_(nullptr),
    incoming_boundary_node_map_(),
    cell_to_outgoing_boundary_nodes_(),
    cell_to_incoming_nonlocal_nodes_(),
    cell_to_outgoing_nonlocal_nodes_()
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

} // namespace opensn