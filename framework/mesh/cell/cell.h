// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/data_types.h"
#include "framework/math/vector3.h"
#include <tuple>
#include <vector>
#include <cstdint>

// Appending cell types to namespace
namespace opensn
{

class MeshContinuum;

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

/// Provides the text name associated with a cell type.
std::string CellTypeName(CellType type);

/**
 * In this paradigm a face is an object which largely is considered to be planar (meaning all the
 * vertices lay in the same plane).
 */
class CellFace
{
public:
  /// Determines the neighbor's partition and whether it's local or not.
  bool IsNeighborLocal(const MeshContinuum* grid) const;

  /// Determines the neighbor's partition.
  int GetNeighborPartitionID(const MeshContinuum* grid) const;

  /// Determines the neighbor's local id.
  uint64_t GetNeighborLocalID(const MeshContinuum* grid) const;

  /// Determines the neighbor's associated face.
  int GetNeighborAdjacentFaceIndex(const MeshContinuum* grid) const;

  /// Computes the geometric info on the face.
  void ComputeGeometricInfo(const MeshContinuum* grid, unsigned int f);

  /// Serializes a face into a vector of bytes.
  ByteArray Serialize() const;

  /// Deserializes a face from a set of raw data
  static CellFace DeSerialize(const ByteArray& raw, size_t& address);

  /// Provides string information of the face.
  std::string ToString() const;

  /// Flag indicating whether face has a neighbor
  bool has_neighbor = false;
  /// If face has neighbor, contains the global_id, otherwise, contains boundary_id.
  uint64_t neighbor_id = 0;

  /// The average/geometric normal
  Vector3 normal;
  /// The face centroid
  Vector3 centroid;
  /// The area of the face
  double area = 0.0;

  /// A list of the vertices
  std::vector<uint64_t> vertex_ids;
};

/// Generic mesh cell object
class Cell
{
public:
  Cell(const Cell& other) = default;
  Cell(Cell&& other) noexcept = default;

  explicit Cell(CellType cell_type, CellType cell_sub_type);

  virtual ~Cell() = default;

  Cell& operator=(const Cell& other);

  CellType GetType() const { return cell_type_; }
  CellType GetSubType() const { return cell_sub_type_; }

  /// Computes the geometric info on the cell.
  void ComputeGeometricInfo(const MeshContinuum* grid);

  /// Serializes a cell into a vector of bytes.
  ByteArray Serialize() const;

  /// Deserializes a cell from a vector of bytes.
  static Cell DeSerialize(const ByteArray& raw, size_t& address);

  /// Provides string information of the cell.
  std::string ToString() const;

  uint64_t global_id = 0;
  uint64_t local_id = 0;
  uint64_t partition_id = 0;
  int material_id = -1;

  Vector3 centroid;
  double volume = 0.0;

  std::vector<uint64_t> vertex_ids;
  std::vector<CellFace> faces;

private:
  /// Primary type, i.e. SLAB, POLYGON, POLYHEDRON
  const CellType cell_type_;
  /// Subtype i.e. SLAB, QUADRILATERAL, HEXAHEDRON
  const CellType cell_sub_type_;
};

} // namespace opensn
