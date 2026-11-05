// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace opensn
{
class FieldFunctionGridBased;

namespace coupling
{

/// Cell-integrated power and geometry for one mesh cell.
struct OwnedCellPower
{
  std::uint64_t global_id = 0;
  unsigned int block_id = 0;
  Vector3 centroid;
  /// Physical three-dimensional volume, including any transverse extent of a 1D or 2D mesh.
  double volume = 0.0;
  /// Cell-integrated power in watts.
  double power = 0.0;
};

/// Cell data passed to a volume-coupling participant.
class CellPowerData
{
public:
  /**
   * Validate and store cell records. `dimension` is the spatial dimension (1, 2, or 3) of
   * the Cartesian mesh the cells came from.
   */
  CellPowerData(unsigned int dimension, std::vector<OwnedCellPower> cells);

  unsigned int GetDimension() const { return dimension_; }
  const std::vector<OwnedCellPower>& GetCells() const { return cells_; }
  double GetLocalPower() const;

  /**
   * Replace cell power density by its global volume-weighted average in each group while
   * preserving every group's integrated power. A group with negative total power is
   * rejected.
   */
  void ConservativelyAverageByGroup(const std::vector<std::int64_t>& local_group_ids);

private:
  unsigned int dimension_;
  std::vector<OwnedCellPower> cells_;
};

/**
 * Integrate one field-function component over each cell and convert to watts.
 *
 * The input field is the energy deposited per unit volume, with volume in the mesh length unit, for
 * the source of the transport input. `source_rate` scales that source to the physical one: physical
 * source particles per second divided by the input's total source rate for the physical body (for
 * 1D and 2D meshes, the input's rate per unit area or depth times `transverse_measure`).
 * `joules_per_energy_unit` converts the field's energy unit to joules, so each cell power is in
 * watts. The mesh must use Cartesian coordinates. A 1D or 2D mesh represents a slab or prism of
 * finite transverse extent: `transverse_measure` is the cross-sectional area of a 1D mesh or the
 * out-of-plane depth of a 2D mesh, in the mesh length unit, and is required for those meshes and
 * rejected for 3D meshes. When `block_ids` is nonempty, only cells in those blocks are included;
 * each listed block must contain at least one cell globally. The returned cell order is the
 * mesh's owned-cell order and contains no ghosts.
 */
CellPowerData MakeCellPowerData(const FieldFunctionGridBased& field_function,
                                std::size_t component,
                                double source_rate,
                                double joules_per_energy_unit,
                                std::optional<double> transverse_measure = std::nullopt,
                                const std::vector<unsigned int>& block_ids = {});

/**
 * Validate cell power for transfer and return its global sum.
 *
 * The global coupling mesh must be nonempty, every rank must report the same dimension, and every
 * cell-integrated power must be nonnegative.
 */
double GetValidatedGlobalPower(const CellPowerData& cell_data);

} // namespace coupling
} // namespace opensn
