#ifndef CELL_MAPPING_FE_PWL_BASE_H
#define CELL_MAPPING_FE_PWL_BASE_H

#include <utility>

#include "mesh/chi_mesh.h"

#include "CellMapping.h"

// ###################################################################
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

  static std::vector<std::vector<int>> MakeFaceNodeMapping(const chi_mesh::Cell& cell);
};

} // namespace chi_math::cell_mapping

#endif