#include "AAH_FLUDSCommonData.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/MeshContinuum/chi_grid_face_histogram.h"

#include "chi_runtime.h"
#include "chi_log.h"

typedef std::vector<std::pair<int, short>> LockBox;

namespace chi_mesh::sweep_management
{

void
AAH_FLUDSCommonData::SlotDynamics(const chi_mesh::Cell& cell,
                                  const SPDS& spds,
                                  const GridFaceHistogram& grid_face_histogram,
                                  std::vector<std::vector<std::pair<int, short>>>& lock_boxes,
                                  std::vector<std::pair<int, short>>& delayed_lock_box,
                                  std::set<int>& location_boundary_dependency_set)
{
  const chi_mesh::MeshContinuum& grid = spds.Grid();

  chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  std::vector<short> inco_face_face_category;
  inco_face_face_category.reserve(cell.faces_.size());
  for (int f = 0; f < cell.faces_.size(); f++)
  {
    const CellFace& face = cell.faces_[f];
    const auto& orientation = spds.CellFaceOrientations()[cell.local_id_][f];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (orientation == FaceOrientation::INCOMING)
    {

      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ LOCAL CELL DEPENDENCE
      if (face.IsNeighborLocal(grid))
      {
        size_t num_face_dofs = face.vertex_ids_.size();
        size_t face_categ = grid_face_histogram.MapFaceHistogramBins(num_face_dofs);

        inco_face_face_category.push_back(static_cast<short>(face_categ));

        LockBox& lock_box = lock_boxes[face_categ];

        //========================================== Check if part of cyclic
        //                                           dependency
        bool is_cyclic = false;
        for (auto cyclic_dependency : spds.GetLocalCyclicDependencies())
        {
          int a = cyclic_dependency.first;
          int b = cyclic_dependency.second;
          int c = cell.local_id_;
          int d = face.GetNeighborLocalID(grid);

          if ((a == c) && (b == d))
          {
            is_cyclic = true;
            inco_face_face_category.back() *= -1;
            inco_face_face_category.back() -= 1;
          }

          if ((a == d) && (b == c))
          {
            is_cyclic = true;
            inco_face_face_category.back() *= -1;
            inco_face_face_category.back() -= 1;
          }
        }
        if (is_cyclic) continue;

        //======================================== Find associated face for
        //                                         dof mapping and lock box
        auto ass_face = (short)face.GetNeighborAssociatedFace(grid);

        // Now find the cell (index,face) pair in the lock box and empty slot
        bool found = false;
        for (auto& lock_box_slot : lock_box)
        {
          if ((lock_box_slot.first == face.neighbor_id_) && (lock_box_slot.second == ass_face))
          {
            lock_box_slot.first = -1;
            lock_box_slot.second = -1;
            found = true;
            break;
          }
        }
        if (!found)
        {
          Chi::log.LogAllError() << "Lock-box location not found in call to "
                                 << "InitializeAlphaElements. Local Cell " << cell.local_id_
                                 << " face " << f << " looking for cell "
                                 << face.GetNeighborLocalID(grid) << " face " << ass_face
                                 << " cat: " << face_categ << " omg=" << spds.Omega().PrintS()
                                 << " lbsize=" << lock_box.size();
          Chi::Exit(EXIT_FAILURE);
        }

      } // if local
      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ BOUNDARY DEPENDENCE
      else if (not face.has_neighbor_)
      {
        const chi_mesh::Vector3& face_norm = face.normal_;

        if (face_norm.Dot(ihat) > 0.999) location_boundary_dependency_set.insert(0);
        else if (face_norm.Dot(ihat) < -0.999)
          location_boundary_dependency_set.insert(1);
        else if (face_norm.Dot(jhat) > 0.999)
          location_boundary_dependency_set.insert(2);
        else if (face_norm.Dot(jhat) < -0.999)
          location_boundary_dependency_set.insert(3);
        else if (face_norm.Dot(khat) > 0.999)
          location_boundary_dependency_set.insert(4);
        else if (face_norm.Dot(khat) < -0.999)
          location_boundary_dependency_set.insert(5);
      }
    } // if incident

  } // for f

  auto raw_inco_face_face_category = new short[inco_face_face_category.size()];
  std::copy(
    inco_face_face_category.begin(), inco_face_face_category.end(), raw_inco_face_face_category);

  so_cell_inco_face_face_category.push_back(raw_inco_face_face_category);

  //=================================================== Loop over faces
  //                OUTGOING                            but process
  //                                                    only outgoing faces
  std::vector<int> outb_face_slot_indices;
  std::vector<short> outb_face_face_category;
  outb_face_slot_indices.reserve(cell.faces_.size());
  outb_face_face_category.reserve(cell.faces_.size());
  for (int f = 0; f < cell.faces_.size(); f++)
  {
    const CellFace& face = cell.faces_[f];
    int cell_g_index = cell.global_id_;
    const auto& orientation = spds.CellFaceOrientations()[cell.local_id_][f];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outgoing face
    if (orientation == FaceOrientation::OUTGOING)
    {
      size_t num_face_dofs = face.vertex_ids_.size();
      size_t face_categ = grid_face_histogram.MapFaceHistogramBins(num_face_dofs);

      outb_face_face_category.push_back(static_cast<short>(face_categ));

      LockBox* temp_lock_box = &lock_boxes[face_categ];

      //========================================== Check if part of cyclic
      //                                           dependency
      if (face.IsNeighborLocal(grid))
      {
        for (auto cyclic_dependency : spds.GetLocalCyclicDependencies())
        {
          int a = cyclic_dependency.first;
          int b = cyclic_dependency.second;
          int c = cell.local_id_;
          int d = face.GetNeighborLocalID(grid);

          if ((a == c) && (b == d))
          {
            temp_lock_box = &delayed_lock_box;
            outb_face_face_category.back() *= -1;
            outb_face_face_category.back() -= 1;
          }

          if ((a == d) && (b == c))
          {
            temp_lock_box = &delayed_lock_box;
            outb_face_face_category.back() *= -1;
            outb_face_face_category.back() -= 1;
          }
        }
      }

      LockBox& lock_box = *temp_lock_box;

      //========================================== Check if this face is
      //                                           the max size
      if (num_face_dofs > largest_face) largest_face = static_cast<int>(num_face_dofs);

      //========================================== Find a open slot
      bool slot_found = false;
      for (int k = 0; k < lock_box.size(); k++)
      {
        if (lock_box[k].first < 0)
        {
          outb_face_slot_indices.push_back(k);
          lock_box[k].first = cell_g_index;
          lock_box[k].second = static_cast<short>(f);
          slot_found = true;
          break;
        }
      }

      //========================================= If an open slot was not found
      //                                          push a new one
      if (!slot_found)
      {
        outb_face_slot_indices.push_back(lock_box.size());
        lock_box.push_back(std::pair<int, short>(cell_g_index, f));
      }

      //========================================== Non-local outgoing
      if (face.has_neighbor_ and (not face.IsNeighborLocal(grid)))
      {
        int locJ = face.GetNeighborPartitionID(grid);
        int deplocI = spds.MapLocJToDeplocI(locJ);
        int face_slot = deplocI_face_dof_count[deplocI];

        deplocI_face_dof_count[deplocI] += face.vertex_ids_.size();

        nonlocal_outb_face_deplocI_slot.emplace_back(deplocI, face_slot);

        // The following function is defined below
        AddFaceViewToDepLocI(deplocI, cell_g_index, face_slot, face);

      } // non-local neighbor
    }   // if outgoing

  } // for f

  auto raw_outb_face_slot_indices = new int[outb_face_slot_indices.size()];
  std::copy(
    outb_face_slot_indices.begin(), outb_face_slot_indices.end(), raw_outb_face_slot_indices);

  so_cell_outb_face_slot_indices.push_back(raw_outb_face_slot_indices);

  auto raw_outb_face_face_category = new short[outb_face_face_category.size()];
  std::copy(
    outb_face_face_category.begin(), outb_face_face_category.end(), raw_outb_face_face_category);

  so_cell_outb_face_face_category.push_back(raw_outb_face_face_category);
}

//###################################################################
/**Given a sweep ordering index, the outgoing face counter,
 * the outgoing face dof, this function computes the location
 * of this position's upwind psi in the local upwind psi vector.*/
void
AAH_FLUDSCommonData::AddFaceViewToDepLocI(int deplocI,
                                          int cell_g_index,
                                          int face_slot,
                                          const chi_mesh::CellFace& face)
{
  //======================================== Check if cell is already there
  bool cell_already_there = false;
  for (auto& cell_view : deplocI_cell_views[deplocI])
  {
    if (cell_view.first == cell_g_index)
    {
      cell_already_there = true;
      cell_view.second.emplace_back(face_slot, face.vertex_ids_);
      break;
    }
  }

  //======================================== If the cell is not there yet
  if (!cell_already_there)
  {
    CompactCellView new_cell_view;
    new_cell_view.first = cell_g_index;
    new_cell_view.second.emplace_back(face_slot, face.vertex_ids_);

    deplocI_cell_views[deplocI].push_back(new_cell_view);
  }
}

} // namespace chi_mesh::sweep_management