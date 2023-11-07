#pragma once

#include <utility>

#include "framework/mesh/mesh.h"

#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"

namespace chi_math::cell_mapping
{
/** Base class for all cell piece-wise linear cell-mappings.
 * \ingroup doc_CellMappings*/
class PieceWiseLinearBaseMapping : public CellMapping
{
protected:
public:
  typedef std::vector<double> VecDbl;
  typedef std::vector<chi_mesh::Vector3> VecVec3;

public:
  /** Constructor. */
  PieceWiseLinearBaseMapping(const chi_mesh::MeshContinuum& grid,
                             const chi_mesh::Cell& cell,
                             size_t num_nodes,
                             std::vector<std::vector<int>> face_node_mappings);

protected:
  static std::vector<chi_mesh::Vector3> GetVertexLocations(const chi_mesh::MeshContinuum& grid,
                                                           const chi_mesh::Cell& cell);

  /** This section just determines a mapping of face dofs
   * to cell dofs. This is pretty simple since we can
   * just loop over each face dof then subsequently
   * loop over cell dofs, if the face dof node index equals
   * the cell dof node index then the mapping is assigned.
   *
   * This mapping is not used by any of the methods in
   *     this class but is used by methods requiring the
   *       surface integrals of the shape functions.
   */
  static std::vector<std::vector<int>> MakeFaceNodeMapping(const chi_mesh::Cell& cell);
};

} // namespace chi_math::cell_mapping
