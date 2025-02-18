// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh.h"
#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include <utility>

namespace opensn
{

/**
 * Base class for all cell piece-wise linear cell-mappings.
 *
 * \ingroup doc_CellMappings
 */
class PieceWiseLinearBaseMapping : public CellMapping
{
protected:
public:
  /// Constructor.
  PieceWiseLinearBaseMapping(const std::shared_ptr<MeshContinuum> grid,
                             const Cell& cell,
                             size_t num_nodes,
                             std::vector<std::vector<int>> face_node_mappings);

protected:
  static std::vector<Vector3> GetVertexLocations(const std::shared_ptr<MeshContinuum> grid,
                                                 const Cell& cell);

  /**
   * This section just determines a mapping of face dofs to cell dofs. This is pretty simple since
   * we can just loop over each face dof then subsequently loop over cell dofs, if the face dof node
   * index equals the cell dof node index then the mapping is assigned.
   *
   * This mapping is not used by any of the methods in  this class but is used by methods requiring
   * the surface integrals of the shape functions.
   */
  static std::vector<std::vector<int>> MakeFaceNodeMapping(const Cell& cell);
};

} // namespace opensn
