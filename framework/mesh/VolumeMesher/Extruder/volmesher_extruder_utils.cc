#include "volmesher_extruder.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "mesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_runtime.h"
#include "chi_log.h"


//###################################################################
/** Creates actual z-levels for the input layer specification.*/
void chi_mesh::VolumeMesherExtruder::SetupLayers(int default_layer_count)
{
  //================================================== Create default layers if no
  //                                                   input layers are provided
  if (input_layers_.empty())
  {
    Chi::log.Log0Warning()
      << "VolumeMesherExtruder: No extrusion layers have been specified. "
      << "A default single layer will be used with height 1.0 and a single "
      << "subdivision.";
    double dz = 1.0/default_layer_count;
    for (int k=0; k<=default_layer_count; k++)
    {
      vertex_layers_.push_back(k * dz);
    }
  }
  else
  {
    double last_z=0.0;
    vertex_layers_.push_back(last_z);

    for (const auto& input_layer : input_layers_)
    {
      double dz = input_layer.height/input_layer.sub_divisions;

      for (int k=0; k<input_layer.sub_divisions; k++)
      {
        last_z += dz;
        vertex_layers_.push_back(last_z);
      }
    }
  }

  Chi::log.Log()
    << "VolumeMesherExtruder: Total number of cell layers is "
    << vertex_layers_.size() - 1;
}

//###################################################################
/**Projects a centroid to an extruded equivalent layer.*/
chi_mesh::Vector3 chi_mesh::VolumeMesherExtruder::
  ProjectCentroidToLevel(const chi_mesh::Vector3 &centroid,
                         const size_t level)
{
  double z_location = 0.5 * (vertex_layers_[level] + vertex_layers_[level + 1]);

  auto centroid_projected = centroid;
  centroid_projected.z = z_location;

  return centroid_projected;
}

//###################################################################
/**Computes a cell's partition id based on a centroid.*/
int chi_mesh::VolumeMesherExtruder::
  GetCellKBAPartitionIDFromCentroid(chi_mesh::Vector3& centroid)
{
  int px = options.partition_x;
  int py = options.partition_y;

  chi_mesh::Cell n_gcell(CellType::GHOST, CellType::GHOST);
  n_gcell.centroid_ = centroid;

  auto xyz_partition_indices = GetCellXYZPartitionID(&n_gcell);

  int nxi = std::get<0>(xyz_partition_indices);
  int nyi = std::get<1>(xyz_partition_indices);
  int nzi = std::get<2>(xyz_partition_indices);

  return nzi*px*py + nyi*px + nxi;
}


//###################################################################
/**Determines if a template cell is in the current partition or a
 * direct neighbor to the current partition.*/
bool chi_mesh::VolumeMesherExtruder::
  HasLocalScope(
    const chi_mesh::Cell& template_cell,
    const chi_mesh::MeshContinuum& template_continuum,
    size_t z_level)
{
  //======================================== Check if the template cell
  //                                         is in the current partition
  {
    auto& centroid = template_cell.centroid_;
    auto projected_centroid = ProjectCentroidToLevel(centroid, z_level);
    int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
    if (pid == Chi::mpi.location_id)
      return true;
  }

  const size_t last_z_level = vertex_layers_.size() - 1;
  const size_t z_level_below = z_level - 1;
  const size_t z_level_above = z_level + 1;

  //======================================== Build z-levels to search
  std::vector<size_t> z_levels_to_search;
  z_levels_to_search.reserve(3);

                                   z_levels_to_search.push_back(z_level);
  if (z_level != 0)                z_levels_to_search.push_back(z_level_below);
  if (z_level != (last_z_level-1)) z_levels_to_search.push_back(z_level_above);

  //======================================== Search template cell's
  //                                         lateral neighbors
  const auto& vertex_subs =
    template_unpartitioned_mesh_->GetVertextCellSubscriptions();
  for (uint64_t vid : template_cell.vertex_ids_)
    for (uint64_t cid : vertex_subs[vid])
    {
      if (cid == template_cell.local_id_) continue;

      auto& candidate_cell =
        template_continuum.local_cells[cid];
      auto& cc_centroid = candidate_cell.centroid_;

      for (size_t z : z_levels_to_search)
      {
        auto projected_centroid = ProjectCentroidToLevel(cc_centroid, z);
        int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
        if (pid == Chi::mpi.location_id)
          return true;
      }
    }//for cid

  //======================================== Search template cell's
  //                                         longitudinal neighbors
  for (size_t z : z_levels_to_search)
  {
    if (z == z_level) continue;

    auto projected_centroid = ProjectCentroidToLevel(template_cell.centroid_, z);
    int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
    if (pid == Chi::mpi.location_id)
      return true;
  }//for z

  return false;
}

//###################################################################
/**Makes an extruded cell from a template cell.*/
std::unique_ptr<chi_mesh::Cell> chi_mesh::VolumeMesherExtruder::
  MakeExtrudedCell(const chi_mesh::Cell &template_cell,
                   const chi_mesh::MeshContinuum &grid,
                   size_t z_level,
                   uint64_t cell_global_id,
                   int partition_id,
                   size_t num_template_cells)
{
  const size_t tc_num_verts = template_cell.vertex_ids_.size();

  //========================================= Determine cell sub-type
  CellType extruded_subtype;
  switch (template_cell.SubType())
  {
    case CellType::TRIANGLE:      extruded_subtype = CellType::WEDGE; break;
    case CellType::QUADRILATERAL: extruded_subtype = CellType::HEXAHEDRON; break;
    default:                      extruded_subtype = CellType::POLYHEDRON;
  }

  //========================================= Create polyhedron
  auto cell = std::make_unique<chi_mesh::Cell>(CellType::POLYHEDRON,
                                               extruded_subtype);
  cell->global_id_    = cell_global_id;
  //cell->local_id set when added to mesh
  cell->partition_id_ = partition_id;
  cell->centroid_ = ProjectCentroidToLevel(template_cell.centroid_, z_level);

  //========================================= Populate cell v-indices
  cell->vertex_ids_.reserve(2 * tc_num_verts);
  for (auto tc_vid : template_cell.vertex_ids_)
    cell->vertex_ids_.push_back(tc_vid + z_level * node_z_index_incr_);

  for (auto tc_vid : template_cell.vertex_ids_)
    cell->vertex_ids_.push_back(tc_vid + (z_level + 1) * node_z_index_incr_);


  //========================================= Create side faces
  for (auto& face : template_cell.faces_)
  {
    chi_mesh::CellFace newFace;

    newFace.vertex_ids_.resize(4, -1);
    newFace.vertex_ids_[0] = face.vertex_ids_[0] + z_level * node_z_index_incr_;
    newFace.vertex_ids_[1] = face.vertex_ids_[1] + z_level * node_z_index_incr_;
    newFace.vertex_ids_[2] = face.vertex_ids_[1] + (z_level + 1) * node_z_index_incr_;
    newFace.vertex_ids_[3] = face.vertex_ids_[0] + (z_level + 1) * node_z_index_incr_;

    //Compute centroid
    const auto& v0 = grid.vertices[newFace.vertex_ids_[0]];
    const auto& v1 = grid.vertices[newFace.vertex_ids_[1]];
    const auto& v2 = grid.vertices[newFace.vertex_ids_[2]];
    const auto& v3 = grid.vertices[newFace.vertex_ids_[3]];

    chi_mesh::Vertex vfc = (v0+v1+v2+v3)/4.0;
    newFace.centroid_ = vfc;

    //Compute normal
    chi_mesh::Vector3 va = v0 - vfc;
    chi_mesh::Vector3 vb = v1 - vfc;

    chi_mesh::Vector3 vn = va.Cross(vb);

    newFace.normal_ = (vn / vn.Norm());

    //Set neighbor
    //The side connections have the same connections as the
    //template cell + the z_level specifiers of the layer.
    if (face.has_neighbor_)
    {
      newFace.neighbor_id_ = face.neighbor_id_ +
                            z_level*num_template_cells;
      newFace.has_neighbor_ = true;
    }
    else
      newFace.neighbor_id_ = face.neighbor_id_;

    cell->faces_.push_back(newFace);
  } //for side faces

  //========================================= Create top and bottom faces
  //=============================== Bottom face
  {
    chi_mesh::CellFace newFace;
    //Vertices
    auto vfc = chi_mesh::Vertex(0.0,0.0,0.0);
    newFace.vertex_ids_.reserve(tc_num_verts);
    for (int tv=(static_cast<int>(tc_num_verts)-1); tv>=0; tv--)
    {
      newFace.vertex_ids_.push_back(template_cell.vertex_ids_[tv]
                                    + z_level * node_z_index_incr_);
      const auto& v = grid.vertices[newFace.vertex_ids_.back()];
      vfc = vfc + v;
    }

    //Compute centroid
    vfc = vfc/static_cast<double>(tc_num_verts);
    newFace.centroid_ = vfc;

    //Compute normal
    auto va = grid.vertices[newFace.vertex_ids_[0]] - vfc;
    auto vb = grid.vertices[newFace.vertex_ids_[1]] - vfc;

    auto vn = va.Cross(vb);
    newFace.normal_ = vn / vn.Norm();

    //Set neighbor
    if (z_level==0)
    {
      newFace.neighbor_id_ = zmin_bndry_id;
      newFace.has_neighbor_ = false;
    }
    else
    {
      newFace.neighbor_id_ = template_cell.local_id_ +
        (z_level-1)*num_template_cells;
      newFace.has_neighbor_ = true;
    }

    cell->faces_.push_back(newFace);
  }

  //=============================== Top face
  {
    chi_mesh::CellFace newFace;
    //Vertices
    auto vfc = chi_mesh::Vertex(0.0,0.0,0.0);
    newFace.vertex_ids_.reserve(tc_num_verts);
    for (auto tc_vid : template_cell.vertex_ids_)
    {
      newFace.vertex_ids_.push_back(tc_vid + (z_level + 1) * node_z_index_incr_);
      const auto& v = grid.vertices[newFace.vertex_ids_.back()];
      vfc = vfc + v;
    }

    //Compute centroid
    vfc = vfc/static_cast<double>(tc_num_verts);
    newFace.centroid_ = vfc;

    //Compute normal
    auto va = grid.vertices[newFace.vertex_ids_[0]] - vfc;
    auto vb = grid.vertices[newFace.vertex_ids_[1]] - vfc;

    auto vn = va.Cross(vb);
    newFace.normal_ = vn / vn.Norm();

    //Set neighbor
    if (z_level==(vertex_layers_.size() - 2))
    {
      newFace.neighbor_id_ = zmax_bndry_id;
      newFace.has_neighbor_ = false;
    }
    else
    {
      newFace.neighbor_id_ = template_cell.local_id_ +
        (z_level+1)*num_template_cells;
      newFace.has_neighbor_ = true;
    }

    cell->faces_.push_back(newFace);
  }

  return cell;
}