// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>
#include <vector>

namespace opensn
{
class MeshContinuum;
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
  void Build(const MeshContinuum& fine_grid, const MeshContinuum& coarse_grif);

  /// Identifier for an invalid face index that means a face maps to nothing.
  static const std::size_t invalid_face_index;

  /// Helper struct for storing the mapping to a coarse cell from a fine cell.
  struct CoarseMapping
  {
    /// Constructor. Sizes fine_faces based on the number of faces within the coarse cell.
    CoarseMapping(const Cell& coarse_cell);
    /// The fine cells contained within a coarse cell.
    std::vector<const Cell*> fine_cells;
    /// The fine cell faces contained within each coarse cell face.
    /// Outer index coarse cell face index (size == # of faces in coarse cell)
    /// Inner index is arbitrary and entries are fine Cell -> fine CellFace index
    std::vector<std::vector<std::pair<const Cell*, std::size_t>>> fine_faces;
  };
  /// Helper struct for storing the mapping from a coarse cell to fine cells.
  struct FineMapping
  {
    /// Constructor. Sizes coarse_faces based on the number of faces within the fine cell.
    FineMapping(const Cell& fine_cell);
    /// The coarse cell that the fine cell is contained within.
    const Cell* coarse_cell;
    /// The coarse CellFace index each fine CellFace is contained within (if any)
    std::vector<std::size_t> coarse_faces;
  };

  /// Get the mapping from the given coarse mesh cell.
  const CoarseMapping& GetCoarseMapping(const Cell& coarse_cell) const;
  /// Get the mapping for the given fine mesh cell.
  const FineMapping& GetFineMapping(const Cell& fine_cell) const;

private:
  /// Mapping for coarse cells to fine cells.
  std::unordered_map<const Cell*, CoarseMapping> coarse_to_fine_;
  /// Mapping for fine cells to a coarse cell.
  std::unordered_map<const Cell*, FineMapping> fine_to_coarse_;
};
} // namespace opensn
