// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include "modules/linear_boltzmann_solvers/lbs_problem/outflow/cell_outflow_view.h"
#include <map>
#include <utility>

namespace opensn
{

class MeshContinuum;

/**
 * Contiguous storage manager for local boundary-face outflow tallies.
 *
 * The bank owns the outflow values and builds one CellOutflowView per local cell. Each view stores
 * offsets into the bank's contiguous storage for the cell's boundary faces.
 */
class OutflowBank
{
public:
  OutflowBank() = default;

  /**
   * Construct outflow storage and cell-local views for a mesh.
   * \param grid Local mesh used to discover boundary faces.
   * \param num_groups Number of energy groups stored per boundary face.
   */
  OutflowBank(const MeshContinuum& grid, unsigned int num_groups);

  /**
   * Transfer cell-local outflow views to the caller.
   * \return Views indexed by local cell id.
   * \note After this operation, the bank no longer owns the view collection, but the returned views
   * still refer to this bank's outflow storage.
   */
  std::vector<CellOutflowView> GetCellOutflowViews() { return std::move(views_); }

  /// Return read/write access to the contiguous outflow values.
  std::vector<double>& GetOutflowData() { return outflows_; }

  /**
   * Return the first group offset for a boundary face.
   * \param cell_local_idx Cell index local to the current process.
   * \param face_idx Face index local to the cell.
   * \return Offset of the face's first group value in the contiguous outflow storage.
   * \throw std::out_of_range If the cell-face pair does not identify a boundary face.
   */
  std::uint64_t GetOffset(std::uint32_t cell_local_idx, std::uint32_t face_idx) const;

  /// Return the number of contiguous outflow values.
  std::size_t GetSize() const { return outflows_.size(); }

private:
  /**
   * Group-wise outflow values for all local boundary faces.
   * Energy groups are the most contiguous index.
   */
  std::vector<double> outflows_;

  /**
   * Map from packed cell-face indices to contiguous outflow offsets.
   * Each offset points to the first energy group for the boundary face.
   */
  std::map<std::uint64_t, std::uint64_t> cellface_map_;

  /// Cell-local views into the contiguous outflow storage.
  std::vector<CellOutflowView> views_;
};

} // namespace opensn
