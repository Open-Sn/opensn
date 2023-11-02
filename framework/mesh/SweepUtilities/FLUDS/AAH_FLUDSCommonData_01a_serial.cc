#include "AAH_FLUDSCommonData.h"

namespace chi_mesh::sweep_management
{

// ###################################################################
/**This cell takes a hierarchy of a cell compact view and
 * serializes it for MPI transmission. This is easy since all
 * the values are integers.*/
void
AAH_FLUDSCommonData::SerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                       std::vector<int>& face_indices,
                                       int num_face_dofs)
{
  const size_t num_cells = cell_views.size();

  //======================== First entry is number of face dofs
  face_indices.push_back(num_face_dofs);

  //======================== Second entry is amount of cells
  face_indices.push_back(static_cast<int>(num_cells));

  //======================== Third entry is negative global cell index
  // Each time a negative entry occurs it denotes a cell face but
  // the actual number is -cell_g_index-1. The offset is necessary
  // for evaluating the negative. The offset is restored during the
  // deserialization process.
  // It is followed by a positive number which is the store location
  // of the face
  for (size_t c = 0; c < num_cells; c++)
  {
    int glob_index = -cell_views[c].first - 1;

    std::vector<CompactFaceView>& cell_face_views = cell_views[c].second;

    size_t num_faces = cell_face_views.size();
    for (size_t f = 0; f < num_faces; f++)
    {
      face_indices.push_back(glob_index);
      face_indices.push_back(cell_face_views[f].first);
      std::vector<uint64_t>& face_vertices = cell_face_views[f].second;

      size_t num_verts = face_vertices.size();
      for (int fi = 0; fi < num_verts; fi++)
      {
        face_indices.push_back(static_cast<int>(face_vertices[fi]));
      }
    }
  }
}

// ###################################################################
/**Deserializes face indices.*/
void
AAH_FLUDSCommonData::DeSerializeCellInfo(std::vector<CompactCellView>& cell_views,
                                         std::vector<int>* face_indices,
                                         int& num_face_dofs)
{
  num_face_dofs = (*face_indices)[0];
  int num_cells = (*face_indices)[1];

  cell_views.resize(num_cells);

  int k = 2;
  int last_cell = -1;
  int c = -1; // cell counter
  int f = -1;
  int v = -1;
  while (k < face_indices->size())
  {
    int entry = (*face_indices)[k];
    //================================= Cell/Face indicator
    if (entry < 0)
    {
      if (-entry != last_cell)
      {
        cell_views.emplace_back();
        c++;
        cell_views[c].first = -entry - 1;

        cell_views[c].second.emplace_back();
        f = 0;

        v = 0;
        last_cell = -entry;

        cell_views[c].second[f].first = (*face_indices)[k + 1];
        k++;
      }
      else
      {
        cell_views[c].second.emplace_back();
        f++;
        v = 0;

        cell_views[c].second[f].first = (*face_indices)[k + 1];
        k++;
      }
    }
    //================================= Face vertex
    else
    {
      cell_views[c].second[f].second.push_back(entry);
      v++;
    }
    k++;
  } // while k
}

} // namespace chi_mesh::sweep_management