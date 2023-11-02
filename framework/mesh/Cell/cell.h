#ifndef CHI_CELL_H
#define CHI_CELL_H

#include "../chi_mesh.h"
#include "data_types/chi_data_types.h"
#include <tuple>

// Appending cell types to namespace
namespace chi_mesh
{
enum class CellType
{
  GHOST = 0,
  SLAB = 1,

  TRIANGLE = 4,
  QUADRILATERAL = 5,
  POLYGON = 6,

  TETRAHEDRON = 7,
  HEXAHEDRON = 8,
  WEDGE = 9,
  PYRAMID = 10,
  POLYHEDRON = 20,

  POINT = 99
};

std::string CellTypeName(CellType type);

//######################################################### Class def
/** In this paradigm a face is an object which largely
 * is considered to be planar (meaning all the vertices
 * lay in the same plane).*/
class CellFace
{
public:
  std::vector<uint64_t> vertex_ids_; /// A list of the vertices
  Normal normal_;                    ///< The average/geometric normal
  Vertex centroid_;                  ///< The face centroid
  bool has_neighbor_ = false;        ///< Flag indicating whether face has a neighbor
  uint64_t neighbor_id_ = 0;         ///< If face has neighbor, contains the global_id.
                                     ///< Otherwise contains boundary_id.

public:
  bool IsNeighborLocal(const chi_mesh::MeshContinuum& grid) const;
  int GetNeighborPartitionID(const chi_mesh::MeshContinuum& grid) const;
  uint64_t GetNeighborLocalID(const chi_mesh::MeshContinuum& grid) const;
  int GetNeighborAssociatedFace(const chi_mesh::MeshContinuum& grid) const;

public:
  double ComputeFaceArea(const chi_mesh::MeshContinuum& grid) const;

  chi_data_types::ByteArray Serialize() const;
  static CellFace DeSerialize(const chi_data_types::ByteArray& raw, size_t& address);
  std::string ToString() const;

  void RecomputeCentroid(const chi_mesh::MeshContinuum& grid);
};

//######################################################### Class def
/**Generic mesh cell object*/
class Cell
{
private:
  const CellType cell_type_;     ///< Primary type, i.e. SLAB, POLYGON, POLYHEDRON
  const CellType cell_sub_type_; ///< Sub-type i.e. SLAB, QUADRILATERAL, HEXAHEDRON

public:
  uint64_t global_id_ = 0;
  uint64_t local_id_ = 0;
  uint64_t partition_id_ = 0;
  Vertex centroid_;
  int material_id_ = -1;

  std::vector<uint64_t> vertex_ids_;
  std::vector<CellFace> faces_;

public:
  Cell(const Cell& other);
  Cell(Cell&& other) noexcept;
  explicit Cell(CellType in_cell_type, CellType in_cell_sub_type)
    : cell_type_(in_cell_type), cell_sub_type_(in_cell_sub_type)
  {
  }

  Cell& operator=(const Cell& other);

  virtual ~Cell() = default;

public:
  CellType Type() const { return cell_type_; }
  CellType SubType() const { return cell_sub_type_; }

  chi_data_types::ByteArray Serialize() const;
  static Cell DeSerialize(const chi_data_types::ByteArray& raw, size_t& address);
  std::string ToString() const;

  void RecomputeCentroidsAndNormals(const chi_mesh::MeshContinuum& grid);
};

} // namespace chi_mesh

#endif