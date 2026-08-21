// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>
#include <vector>
#include <memory>

namespace opensn
{
class Mesh;
class Cell;

/**
 * Produces a mapping between a fine mesh and a coarse mesh.
 * Each fine mesh cell face must be fully contained within a
 * coarse mesh face. The meshes can be the same.
 */
class MeshMapping
{
public:
  MeshMapping() = default;

  /// Builds the mapping.
  void Build(const std::shared_ptr<Mesh>& fine_grid, const std::shared_ptr<Mesh>& coarse_grid);

  /// Helper struct for storing the mapping to a coarse cell from a fine cell.
  struct CoarseMapping
  {
    /// Constructor. Sizes fine_faces based on the number of faces within the coarse cell.
    explicit CoarseMapping(size_t num_faces);

    /// The fine cells contained within a coarse cell.
    std::vector<std::uint32_t> fine_cell_local_ids;
    /// The fine cell faces contained within each coarse cell face.
    /// Outer index coarse cell face index (size == # of faces in coarse cell)
    /// Inner index is arbitrary and entries are fine Cell -> fine CellFace index
    std::vector<std::vector<std::pair<const Cell*, size_t>>> fine_faces;
  };

  /// Helper struct for storing the mapping from a coarse cell to fine cells.
  struct FineMapping
  {
    /// Constructor. Sizes coarse_faces based on the number of faces within the fine cell.
    explicit FineMapping(size_t num_faces);

    /// The coarse cell that the fine cell is contained within.
    std::uint32_t coarse_cell_local_id;
    /// The coarse CellFace index each fine CellFace is contained within (if any)
    std::vector<size_t> coarse_faces;
  };

  /// Get the mapping from the given coarse mesh cell.
  const CoarseMapping& GetCoarseMapping(std::uint32_t coarse_cell_local_id) const;
  /// Get the mapping for the given fine mesh cell.
  const FineMapping& GetFineMapping(std::uint32_t cell_local_id) const;

private:
  /// Mapping for coarse cells to fine cells.
  std::unordered_map<std::uint32_t, CoarseMapping> coarse_to_fine_;
  /// Mapping for fine cells to a coarse cell : local cell id of a fine cell -> cell mapping
  std::unordered_map<std::uint32_t, FineMapping> fine_to_coarse_;

public:
  /// Identifier for an invalid face index that means a face maps to nothing.
  static const size_t invalid_face_index;
};
} // namespace opensn
