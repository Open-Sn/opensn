// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh.h"
#include "framework/data_types/data_types.h"
#include <tuple>

// Appending cell types to namespace
namespace opensn
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

/**Provides the text name associated with a cell type.*/
std::string CellTypeName(CellType type);

/** In this paradigm a face is an object which largely
 * is considered to be planar (meaning all the vertices
 * lay in the same plane).*/
class CellFace
{
public:
  std::vector<uint64_t> vertex_ids_; /// A list of the vertices
  Vector3 normal_;                   ///< The average/geometric normal
  Vector3 centroid_;                 ///< The face centroid
  bool has_neighbor_ = false;        ///< Flag indicating whether face has a neighbor
  uint64_t neighbor_id_ = 0;         ///< If face has neighbor, contains the global_id.
                                     ///< Otherwise contains boundary_id.

public:
  /**Determines the neighbor's partition and whether its local or not.*/
  bool IsNeighborLocal(const MeshContinuum& grid) const;
  /**Determines the neighbor's partition.*/
  int GetNeighborPartitionID(const MeshContinuum& grid) const;
  /**Determines the neighbor's local id.*/
  uint64_t GetNeighborLocalID(const MeshContinuum& grid) const;
  /**Determines the neighbor's associated face.*/
  int GetNeighborAssociatedFace(const MeshContinuum& grid) const;

public:
  /**Computes the face area.*/
  double ComputeFaceArea(const MeshContinuum& grid) const;

  /**Serializes a face into a vector of bytes.*/
  ByteArray Serialize() const;
  /**Deserializes a face from a set of raw data*/
  static CellFace DeSerialize(const ByteArray& raw, size_t& address);
  /**Provides string information of the face.*/
  std::string ToString() const;

  /**Recomputes the face centroid assuming the mesh vertices have been transformed.*/
  void RecomputeCentroid(const MeshContinuum& grid);
};

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
  Vector3 centroid_;
  int material_id_ = -1;

  std::vector<uint64_t> vertex_ids_;
  std::vector<CellFace> faces_;

public:
  /**Copy constructor*/
  Cell(const Cell& other);
  /**Move constructor*/
  Cell(Cell&& other) noexcept;
  explicit Cell(CellType cell_type, CellType cell_sub_type)
    : cell_type_(cell_type), cell_sub_type_(cell_sub_type)
  {
  }

  /**Copy operator.*/
  Cell& operator=(const Cell& other);

  virtual ~Cell() = default;

public:
  CellType Type() const { return cell_type_; }
  CellType SubType() const { return cell_sub_type_; }

  /**Serializes a cell into a vector of bytes.*/
  ByteArray Serialize() const;
  /**Deserializes a cell from a vector of bytes.*/
  static Cell DeSerialize(const ByteArray& raw, size_t& address);
  /**Provides string information of the cell.*/
  std::string ToString() const;

  /**Recomputes the cell centroid and all face centroids assuming
   * the mesh vertices have been transformed.*/
  void RecomputeCentroidsAndNormals(const MeshContinuum& grid);
};

} // namespace opensn
