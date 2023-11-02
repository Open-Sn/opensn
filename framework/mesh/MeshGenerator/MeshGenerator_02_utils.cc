#include "MeshGenerator.h"

namespace chi_mesh
{

// ###################################################################
/**Broadcasts PIDs to other locations.*/
void
MeshGenerator::BroadcastPIDs(std::vector<int64_t>& cell_pids, int root, MPI_Comm communicator)
{
  size_t data_count = Chi::mpi.location_id == root ? cell_pids.size() : 0;

  //======================================== Broadcast data_count to all
  //                                         locations
  MPI_Bcast(&data_count,   // buffer [IN/OUT]
            1,             // count
            MPI_UINT64_T,  // data type
            root,          // root
            communicator); // communicator

  if (Chi::mpi.location_id != root) cell_pids.assign(data_count, 0);

  //======================================== Broadcast partitioning to all
  //                                         locations
  MPI_Bcast(cell_pids.data(),             // buffer [IN/OUT]
            static_cast<int>(data_count), // count
            MPI_LONG_LONG_INT,            // data type
            root,                         // root
            communicator);                // communicator
}

// ###################################################################
/**Determines if a cells needs to be included as a ghost or as a local cell.*/
bool
MeshGenerator::CellHasLocalScope(int location_id,
                                 const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
                                 uint64_t cell_global_id,
                                 const std::vector<std::set<uint64_t>>& vertex_subscriptions,
                                 const std::vector<int64_t>& cell_partition_ids) const
{
  if (replicated_) return true;
  // First determine if the cell is a local cell
  int cell_pid = static_cast<int>(cell_partition_ids[cell_global_id]);
  if (cell_pid == location_id) return true;

  // Now determine if the cell is a ghost cell
  for (uint64_t vid : lwcell.vertex_ids)
    for (uint64_t cid : vertex_subscriptions[vid])
    {
      if (cid == cell_global_id) continue;
      int adj_pid = static_cast<int>(cell_partition_ids[cid]);
      if (adj_pid == location_id) return true;
    }

  return false;
}

// ###################################################################
/**Converts a light-weight cell to a real cell.*/
std::unique_ptr<chi_mesh::Cell>
MeshGenerator::SetupCell(const UnpartitionedMesh::LightWeightCell& raw_cell,
                         uint64_t global_id,
                         uint64_t partition_id,
                         const VertexListHelper& vertices)
{
  auto cell = std::make_unique<chi_mesh::Cell>(raw_cell.type, raw_cell.sub_type);
  cell->centroid_ = raw_cell.centroid;
  cell->global_id_ = global_id;
  cell->partition_id_ = partition_id;
  cell->material_id_ = raw_cell.material_id;

  cell->vertex_ids_ = raw_cell.vertex_ids;

  size_t face_counter = 0;
  for (auto& raw_face : raw_cell.faces)
  {
    chi_mesh::CellFace newFace;

    newFace.has_neighbor_ = raw_face.has_neighbor;
    newFace.neighbor_id_ = raw_face.neighbor;

    newFace.vertex_ids_ = raw_face.vertex_ids;
    auto vfc = chi_mesh::Vertex(0.0, 0.0, 0.0);
    for (auto fvid : newFace.vertex_ids_)
      vfc = vfc + vertices.at(fvid);
    newFace.centroid_ = vfc / static_cast<double>(newFace.vertex_ids_.size());

    if (cell->Type() == CellType::SLAB)
    {
      // A slab face is very easy. If it is the first face
      // the normal is -khat. If it is the second face then
      // it is +khat.
      if (face_counter == 0) newFace.normal_ = chi_mesh::Vector3(0.0, 0.0, -1.0);
      else
        newFace.normal_ = chi_mesh::Vector3(0.0, 0.0, 1.0);
    }
    else if (cell->Type() == CellType::POLYGON)
    {
      // A polygon face is just a line so we can just grab
      // the first vertex and form a vector with the face
      // centroid. The normal is then just khat
      // cross-product with this vector.
      uint64_t fvid = newFace.vertex_ids_[0];
      auto vec_vvc = vertices.at(fvid) - newFace.centroid_;

      newFace.normal_ = chi_mesh::Vector3(0.0, 0.0, 1.0).Cross(vec_vvc);
      newFace.normal_.Normalize();
    }
    else if (cell->Type() == CellType::POLYHEDRON)
    {
      // A face of a polyhedron can itself be a polygon
      // which can be multifaceted. Here we need the
      // average normal over all the facets computed
      // using an area-weighted average.
      const size_t num_face_verts = newFace.vertex_ids_.size();
      double total_area = 0.0;
      for (size_t fv = 0; fv < num_face_verts; ++fv)
      {
        size_t fvp1 = (fv < (num_face_verts - 1)) ? fv + 1 : 0;

        uint64_t fvid_m = newFace.vertex_ids_[fv];
        uint64_t fvid_p = newFace.vertex_ids_[fvp1];

        auto leg_m = vertices.at(fvid_m) - newFace.centroid_;
        auto leg_p = vertices.at(fvid_p) - newFace.centroid_;

        auto vn = leg_m.Cross(leg_p);

        double area = 0.5 * vn.Norm();
        total_area += area;

        newFace.normal_ = newFace.normal_ + area * vn.Normalized();
      }
      newFace.normal_ = newFace.normal_ / total_area;
      newFace.normal_.Normalize();
    }
    ++face_counter;

    cell->faces_.push_back(newFace);
  }

  return cell;
}

} // namespace chi_mesh