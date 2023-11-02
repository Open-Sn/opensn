#include "chi_unpartitioned_mesh.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

#include "chi_mpi.h"

//###################################################################
/**Establishes neighbor connectivity for the light-weight mesh.*/
void chi_mesh::UnpartitionedMesh::BuildMeshConnectivity()
{
  const size_t num_raw_cells = raw_cells_.size();
  const size_t num_raw_vertices = vertices_.size();

  //======================================== Reset all cell neighbors
  int num_bndry_faces = 0;
  for (auto& cell : raw_cells_)
    for (auto& face : cell->faces)
      if (not face.has_neighbor) ++num_bndry_faces;

  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                              << " Number of unconnected faces "
                                 "before connectivity: " << num_bndry_faces;

  Chi::log.Log() << Chi::program_timer.GetTimeString()
                << " Establishing cell connectivity.";

  //======================================== Establish internal connectivity
  // Populate vertex subscriptions to internal cells
  vertex_cell_subscriptions_.resize(num_raw_vertices);
  {
    uint64_t cur_cell_id=0;
    for (const auto& cell : raw_cells_)
    {
      for (auto vid : cell->vertex_ids)
        vertex_cell_subscriptions_.at(vid).insert(cur_cell_id);
      ++cur_cell_id;
    }
  }

  Chi::log.Log() << Chi::program_timer.GetTimeString()
                << " Vertex cell subscriptions complete.";

  // Process raw cells
  {
    uint64_t aux_counter = 0;
    uint64_t cur_cell_id=0;
    for (auto& cell : raw_cells_)
    {
      for (auto& cur_cell_face : cell->faces)
      {
        if (cur_cell_face.has_neighbor) {continue;}
        const std::set<uint64_t> cfvids(cur_cell_face.vertex_ids.begin(),
                                        cur_cell_face.vertex_ids.end());

        std::set<size_t> cells_to_search;
        for (uint64_t vid : cfvids)
          for (uint64_t cell_id : vertex_cell_subscriptions_.at(vid))
            if (cell_id != cur_cell_id)
              cells_to_search.insert(cell_id);

        for (uint64_t adj_cell_id : cells_to_search)
        {
          auto adj_cell = raw_cells_.at(adj_cell_id);

          for (auto& adj_cell_face : adj_cell->faces)
          {
            if (adj_cell_face.has_neighbor) {continue;}
            const std::set<uint64_t> afvids(adj_cell_face.vertex_ids.begin(),
                                            adj_cell_face.vertex_ids.end());

            if (cfvids == afvids)
            {
              cur_cell_face.neighbor = adj_cell_id;
              adj_cell_face.neighbor = cur_cell_id;

              cur_cell_face.has_neighbor = true;
              adj_cell_face.has_neighbor = true;

              goto face_neighbor_found;
            }
          }//for adjacent cell face
        }
        face_neighbor_found:;
      }//for face

      ++cur_cell_id;
      const double fraction_complete = static_cast<double>(cur_cell_id)/
                                       static_cast<double>(num_raw_cells);
      if (fraction_complete >= static_cast<double>(aux_counter+1)*0.1)
      {
        Chi::log.Log() << Chi::program_timer.GetTimeString()
                      << " Surpassing cell " << cur_cell_id
                      << " of " << num_raw_cells
                      << " (" << (aux_counter+1)*10 << "%)";
        ++aux_counter;
      }
    }//for cell
  }

  Chi::log.Log() << Chi::program_timer.GetTimeString()
                << " Establishing cell boundary connectivity.";

  //======================================== Establish boundary connectivity
  // Make list of internal cells on the boundary
  std::vector<LightWeightCell*> internal_cells_on_boundary;
  for (auto& cell : raw_cells_)
  {
    bool cell_on_boundary = false;
    for (auto& face : cell->faces)
      if (not face.has_neighbor)
      { cell_on_boundary = true; break; }

    if (cell_on_boundary) internal_cells_on_boundary.push_back(cell);
  }

  // Populate vertex subscriptions to boundary cells
  std::vector<std::set<uint64_t>> vertex_bndry_cell_subscriptions(vertices_.size());
  {
    uint64_t cur_cell_id=0;
    for (auto& cell : raw_boundary_cells_)
    {
      for (auto vid : cell->vertex_ids)
        vertex_bndry_cell_subscriptions.at(vid).insert(cur_cell_id);
      ++cur_cell_id;
    }
  }

  // Process boundary cells
  for (auto& cell : internal_cells_on_boundary)
    for (auto& face : cell->faces)
    {
      if (face.has_neighbor) continue;
      std::set<uint64_t> cfvids(face.vertex_ids.begin(),
                                face.vertex_ids.end());

      std::set<size_t> cells_to_search;
      for (uint64_t vid : face.vertex_ids)
        for (uint64_t cell_id : vertex_bndry_cell_subscriptions[vid])
          cells_to_search.insert(cell_id);

      for (uint64_t adj_cell_id : cells_to_search)
      {
        auto& adj_cell = raw_boundary_cells_[adj_cell_id];

        std::set<uint64_t> afvids(adj_cell->vertex_ids.begin(),
                                  adj_cell->vertex_ids.end());

        if (cfvids == afvids)
        {
          face.neighbor = adj_cell->material_id;
          break;
        }
      }//for adj_cell_id
    }//for face

  num_bndry_faces = 0;
  for (auto cell : raw_cells_)
    for (auto& face : cell->faces)
      if (not face.has_neighbor) ++num_bndry_faces;

  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                              << " Number of boundary faces "
                                 "after connectivity: " << num_bndry_faces;

  Chi::log.Log() << Chi::program_timer.GetTimeString()
                << " Done establishing cell connectivity.";

}