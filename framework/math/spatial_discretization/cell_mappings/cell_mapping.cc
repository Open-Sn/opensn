// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/data_types/matrix3x3.h"
#include <utility>

namespace opensn
{

CellMapping::CellMapping(const std::shared_ptr<MeshContinuum> grid,
                         const Cell& cell,
                         size_t num_nodes,
                         std::vector<Vector3> node_locations,
                         std::vector<std::vector<int>> face_node_mappings)
  : grid_(grid),
    cell_(cell),
    num_nodes_(num_nodes),
    node_locations_(std::move(node_locations)),
    face_node_mappings_(std::move(face_node_mappings))
{
}

} // namespace opensn
