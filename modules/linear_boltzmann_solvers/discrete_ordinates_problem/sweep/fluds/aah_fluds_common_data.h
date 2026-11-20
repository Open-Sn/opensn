// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <set>
#include <unordered_map>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>

namespace opensn
{

class GridFaceHistogram;
class Cell;
class CellFace;

// face_slot index, vertex ids
using CompactFaceView = std::pair<int64_t, std::vector<uint64_t>>;

// cell_global_id, faces
using CompactCellView = std::pair<uint64_t, std::vector<CompactFaceView>>;

class SPDS;

class AAH_FLUDSCommonData : public FLUDSCommonData
{
public:
  /// Construct alpha phase only (local work, no MPI). Call FinalizeBeta() afterwards.
  /// MakeAlpha() prevents us from holding an unfinalized object.
  static std::unique_ptr<AAH_FLUDSCommonData>
  MakeAlpha(const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
            const SPDS& spds,
            const GridFaceHistogram& grid_face_histogram);

private:
  explicit AAH_FLUDSCommonData(const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
                               const SPDS& spds,
                               const GridFaceHistogram& grid_face_histogram);

public:
  /// Complete SPDS construction. Exchange inter-rank face data via MPI. This routine must be
  /// called sequentially across all SPDS in the same order on all ranks after the constructor.
  void FinalizeBeta(const SPDS& spds);

  bool IsFinalized() const { return finalized_; }

protected:
  friend class AAH_FLUDS;
  size_t largest_face_ = 0;
  /// Number of face categories
  size_t num_face_categories_ = 0;
  /// Group-angle-faceDOF stride per cat
  std::vector<size_t> local_psi_stride_;
  /// Number of faces in each cat
  std::vector<size_t> local_psi_max_elements_;
  /// Group-angle-faceDOF stride delayed cat
  size_t delayed_local_psi_stride_ = 0;
  /// Number of faces in delayed cat
  size_t delayed_local_psi_max_elements_ = 0;

  /**
   * local_psi_n_block_stride[fc]. Given face category fc, the value is
   * total number of faces that store information in this category's buffer
   * per angle
   */
  std::vector<size_t> local_psi_n_block_stride_;
  std::vector<size_t> local_psi_Gn_block_strideG_;
  size_t delayed_local_psi_Gn_block_stride_ = 0;
  size_t delayed_local_psi_Gn_block_strideG_ = 0;

  /**
   * This is a small vector [deplocI] that holds the number of face dofs for each dependent
   * location.
   */
  std::vector<uint64_t> deplocI_face_dof_count_;

  /**
   * This is a vector [dependent_location][unordered_cell_index]
   * that holds an AlphaPair. AlphaPair-first is the cell's global_id
   * and AlphaPair-second holds a number of BetaPairs. Each BetaPair
   * comprises BetaPair-first = face_slot_index (the location of this
   * faces data in the psi vector, and then a vector of vertex indexes
   * that can be used for dof_mapping.
   * Filled during slot-dynamics.
   * Cleared after beta-pass.
   */
  std::vector<std::vector<CompactCellView>> deplocI_cell_views_;
  std::vector<std::unordered_map<uint64_t, size_t>> deploc_i_cell_idx_;
  bool finalized_ = false;

  /**
   * The outgoing faces of the cell with sweep order index s are entries outb_face_offsets_[s] to
   * outb_face_offsets_[s + 1] - 1 of the outgoing-face arrays, in the order of the cell's faces.
   * The same holds for the incoming faces with a local neighbor and inco_face_offsets_. Every
   * direction has its own sweep ordering and so its own copy of these arrays, so they are stored
   * flat rather than as a vector per cell.
   */
  std::vector<std::uint32_t> outb_face_offsets_;

  /// Lock-box slot of each outgoing face, where its first face node stores its data.
  std::vector<std::uint32_t> outb_face_slots_;

  /**
   * Face category of each outgoing face, that is, the local psi vector that holds faces of its
   * size. A delayed face, whose downwind cell reads lagged data, stores -category - 1.
   */
  std::vector<short> outb_face_categories_;

  /// Offsets of the incoming faces with a local neighbor, as outb_face_offsets_.
  std::vector<std::uint32_t> inco_face_offsets_;

  /// Face category of each incoming face, encoded as outb_face_categories_.
  std::vector<short> inco_face_categories_;

  /// Lock-box slot where the upwind cell stores the data of each incoming face.
  std::vector<std::uint32_t> inco_face_slots_;

  /**
   * The grid nodal mapping of each incoming face, which maps its nodes to the nodes of the upwind
   * face. The pointers address the face node mappings of the problem's grid nodal mappings, which
   * are built once with the problem and outlive the sweep data built from them.
   */
  std::vector<const short*> inco_face_upwind_dofs_;

private:
  void InitializeAlphaElements(const SPDS& spds, const GridFaceHistogram& grid_face_histogram);

  void SlotDynamics(std::uint32_t cell_local_id,
                    const SPDS& spds,
                    const GridFaceHistogram& grid_face_histogram,
                    std::vector<std::vector<std::pair<std::optional<uint64_t>, short>>>& lock_boxes,
                    std::vector<std::pair<std::optional<uint64_t>, short>>& delayed_lock_box);

  /**
   * Given a sweep ordering index, the outgoing face counter, the outgoing face dof, this function
   * computes the location of this position's upwind psi in the local upwind psi vector.
   */
  void AddFaceViewToDepLocI(int deplocI,
                            uint64_t cell_g_index,
                            uint64_t face_slot,
                            const std::vector<uint64_t>& face_vertex_ids);

  void LocalIncidentMapping(std::uint32_t cell_local_id,
                            const SPDS& spds,
                            std::vector<uint64_t>& local_so_cell_mapping);

  void InitializeBetaElements(const SPDS& spds, int tag_index = 0);

  void
  NonLocalIncidentMapping(std::uint32_t cell_local_id,
                          const SPDS& spds,
                          const std::vector<std::unordered_map<uint64_t, size_t>>& preloc_i_idx,
                          const std::vector<std::unordered_map<uint64_t, size_t>>& dpreloc_i_idx);

  /**
   * This is a vector [non_local_outgoing_face_count]  that maps a face to a dependent location and
   * associated slot index
   */
  std::vector<std::pair<int, int>> nonlocal_outb_face_deplocI_slot_;

  /**
   * This is a vector [predecessor_location][unordered_cell_index]
   * that holds an AlphaPair. AlphaPair-first is the cell's global_id
   * and AlphaPair-second holds a number of BetaPairs. Each BetaPair
   * comprises BetaPair-first = face_slot_index (the location of this
   * faces data in the psi vector, and then a vector of vertex indexes
   * that can be used for dof_mapping.
   * Filled in beta-pass
   * Cleared after beta-pass.
   */
  std::vector<std::vector<CompactCellView>> prelocI_cell_views_;
  std::vector<std::vector<CompactCellView>> delayed_prelocI_cell_views_;

  /**
   * This is a small vector [prelocI] that holds the number of  face dofs for each predecessor
   * location.
   */
  std::vector<int64_t> prelocI_face_dof_count_;
  std::vector<int64_t> delayed_prelocI_face_dof_count_;

  /**
   * This is a vector [nonlocal_inc_face_counter] containing AlphaPairs. AlphaPair-first is the
   * prelocI index and AlphaPair-second is a BetaPair. The BetaPair-first is the slot where
   * the face storage begins and BetaPair-second is a dof mapping
   */
  std::vector<std::pair<int, std::pair<int64_t, std::vector<int64_t>>>>
    nonlocal_inc_face_prelocI_slot_dof_;

  std::vector<std::pair<int, std::pair<int64_t, std::vector<int64_t>>>>
    delayed_nonlocal_inc_face_prelocI_slot_dof_;

public:
  /**
   * This cell takes a hierarchy of a cell compact view and serializes it for MPI transmission.
   * This is easy since all the values are integers.
   */
  static void SerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                std::vector<int64_t>& face_indices,
                                uint64_t num_face_dofs);

  /// Deserializes face indices.
  static void DeSerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                  std::vector<int64_t>* face_indices,
                                  int64_t& num_face_dofs);
};

} // namespace opensn
