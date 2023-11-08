#pragma once

#include "framework/mesh/SweepUtilities/FLUDS/FLUDSCommonData.h"
#include <cstddef>
#include <cstdint>
#include <set>

namespace chi_mesh
{
class GridFaceHistogram;
class Cell;
class CellFace;
} // namespace chi_mesh

namespace chi_mesh::sweep_management
{

// face_slot index, vertex ids
typedef std::pair<int, std::vector<uint64_t>> CompactFaceView;

// cell_global_id, faces
typedef std::pair<int, std::vector<CompactFaceView>> CompactCellView;

class SPDS;

class AAH_FLUDSCommonData : public FLUDSCommonData
{
public:
  struct INCOMING_FACE_INFO
  {
    int slot_address = 0;
    short* upwind_dof_mapping = nullptr;

    void Setup(const std::pair<int, std::vector<short>>& input)
    {
      slot_address = input.first;
      upwind_dof_mapping = new short[input.second.size()];
      std::copy(input.second.begin(), input.second.end(), upwind_dof_mapping);
    }

    ~INCOMING_FACE_INFO() { delete[] upwind_dof_mapping; }
  }; // TODO: Make common
public:
  explicit AAH_FLUDSCommonData(const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
                               const SPDS& spds,
                               const chi_mesh::GridFaceHistogram& grid_face_histogram);

protected:
  friend class AAH_FLUDS;
  int largest_face = 0;
  /// Number of face categories
  size_t num_face_categories = 0;
  /// Group-angle-faceDOF stride per cat
  std::vector<size_t> local_psi_stride;
  /// Number of faces in each cat
  std::vector<size_t> local_psi_max_elements;
  /// Group-angle-faceDOF stride delayed cat
  size_t delayed_local_psi_stride = 0;
  /// Number of faces in delayed cat
  size_t delayed_local_psi_max_elements = 0;

  /// local_psi_n_block_stride[fc]. Given face category fc, the value is
  /// total number of faces that store information in this category's buffer
  /// per angle
  std::vector<size_t> local_psi_n_block_stride;
  std::vector<size_t> local_psi_Gn_block_strideG;
  size_t delayed_local_psi_Gn_block_stride = 0;
  size_t delayed_local_psi_Gn_block_strideG = 0;

  /// Very small vector listing the boundaries this location depends on
  std::vector<int> boundary_dependencies;

  /// This is a small vector [deplocI] that holds the number of
  /// face dofs for each dependent location.
  std::vector<int> deplocI_face_dof_count;

  /// This is a vector [dependent_location][unordered_cell_index]
  /// that holds an AlphaPair. AlphaPair-first is the cell's global_id
  /// and AlphaPair-second holds a number of BetaPairs. Each BetaPair
  /// comprises BetaPair-first = face_slot_index (the location of this
  /// faces data in the psi vector, and then a vector of vertex indexes
  /// that can be used for dof_mapping.
  /// Filled during slot-dynamics.
  /// Cleared after beta-pass.
  std::vector<std::vector<CompactCellView>> deplocI_cell_views;

  /// This is a vector [cell_sweep_order_index][outgoing_face_count]
  /// which holds the slot address in the local psi vector where the first
  /// face dof will store its data
  std::vector<int*> so_cell_outb_face_slot_indices;

  /// This is a vector [cell_sweep_order_index][outgoing_face_count]
  /// which holds the face categorization for the face. i.e. the local
  /// psi vector that hold faces of the same category.
  std::vector<short*> so_cell_outb_face_face_category;

  /// This is a vector [cell_sweep_order_index][incoming_face_count]
  /// which holds the face categorization for the face. i.e. the local
  /// psi vector that hold faces of the same category.
  std::vector<short*> so_cell_inco_face_face_category;

  /// This is a vector [cell_sweep_order_index][incoming_face_count]
  /// that will hold a structure. struct.slot_address holds the slot address
  /// where this face's upwind data is stored. struct.upwind_dof_mapping is
  /// a mapping of each of this face's dofs to the upwinded face's dofs
private:
  std::vector<INCOMING_FACE_INFO*> so_cell_inco_face_dof_indices;

  /// This is a vector [non_local_outgoing_face_count]
  /// that maps a face to a dependent location and associated slot index
  std::vector<std::pair<int, int>> nonlocal_outb_face_deplocI_slot;

private:
  /// This is a vector [predecessor_location][unordered_cell_index]
  /// that holds an AlphaPair. AlphaPair-first is the cell's global_id
  /// and AlphaPair-second holds a number of BetaPairs. Each BetaPair
  /// comprises BetaPair-first = face_slot_index (the location of this
  /// faces data in the psi vector, and then a vector of vertex indexes
  /// that can be used for dof_mapping.
  /// Filled in beta-pass
  /// Cleared after beta-pass.
  std::vector<std::vector<CompactCellView>> prelocI_cell_views;
  std::vector<std::vector<CompactCellView>> delayed_prelocI_cell_views;

  /// This is a small vector [prelocI] that holds the number of
  /// face dofs for each predecessor location.
  std::vector<int> prelocI_face_dof_count;
  std::vector<int> delayed_prelocI_face_dof_count;

  /// This is a vector [nonlocal_inc_face_counter] containing
  /// AlphaPairs. AlphaPair-first is the prelocI index and
  /// AlphaPair-second is a BetaPair. The BetaPair-first is the slot where
  /// the face storage begins and BetaPair-second is a dof mapping
  std::vector<std::pair<int, std::pair<int, std::vector<int>>>> nonlocal_inc_face_prelocI_slot_dof;

  std::vector<std::pair<int, std::pair<int, std::vector<int>>>>
    delayed_nonlocal_inc_face_prelocI_slot_dof;

  void InitializeAlphaElements(const SPDS& spds, const GridFaceHistogram& grid_face_histogram);

  void SlotDynamics(const chi_mesh::Cell& cell,
                    const SPDS& spds,
                    const GridFaceHistogram& grid_face_histogram,
                    std::vector<std::vector<std::pair<int, short>>>& lock_boxes,
                    std::vector<std::pair<int, short>>& delayed_lock_box,
                    std::set<int>& location_boundary_dependency_set);

  /**
   * Given a sweep ordering index, the outgoing face counter, the outgoing face dof, this function
   * computes the location of this position's upwind psi in the local upwind psi vector.
   */
  void AddFaceViewToDepLocI(int deplocI,
                            int cell_g_index,
                            int face_slot,
                            const chi_mesh::CellFace& face);

  void LocalIncidentMapping(const chi_mesh::Cell& cell,
                            const SPDS& spds,
                            std::vector<int>& local_so_cell_mapping);

  void InitializeBetaElements(const SPDS& spds, int tag_index = 0);

  /**
   * This cell takes a hierarchy of a cell compact view and serializes it for MPI transmission.
   * This is easy since all the values are integers.
   */
  static void SerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                std::vector<int>& face_indices,
                                int num_face_dofs);

  /**
   * Deserializes face indices.
   */
  static void DeSerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                  std::vector<int>* face_indices,
                                  int& num_face_dofs);

  void NonLocalIncidentMapping(const chi_mesh::Cell& cell, const SPDS& spds);
};

} // namespace chi_mesh::sweep_management
