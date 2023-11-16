#include "framework/mesh/volume_mesher/extruder/volmesher_extruder.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/surface_mesher/surface_mesher.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/utils/timer.h"
#include "framework/console/console.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi.h"
#include <iostream>

namespace opensn
{

void
VolumeMesherExtruder::CreateLocalNodes(MeshContinuum& template_grid, MeshContinuum& grid)
{
  // For each layer
  std::set<uint64_t> vertex_ids_with_local_scope;
  for (size_t iz = 0; iz < (vertex_layers_.size() - 1); iz++)
  {
    for (const auto& template_cell : template_grid.local_cells)
    {
      // Check template cell type
      if (template_cell.Type() != CellType::POLYGON)
        throw std::logic_error("Extruder::CreateLocalNodes: "
                               "Template cell error. Not of base type POLYGON");

      bool has_local_scope = HasLocalScope(template_cell, template_grid, iz);

      if (has_local_scope)
      {
        auto& vertex_set = vertex_ids_with_local_scope;
        for (auto tc_vid : template_cell.vertex_ids_)
          vertex_set.insert(tc_vid + iz * node_z_index_incr_);

        for (auto tc_vid : template_cell.vertex_ids_)
          vertex_set.insert(tc_vid + (iz + 1) * node_z_index_incr_);
      }
    } // for template cell
  }   // for layer

  // Now add all nodes that are local or neighboring
  uint64_t vid = 0;
  for (auto layer_z_level : vertex_layers_)
  {
    for (auto& id_vertex : template_grid.vertices)
    {
      const auto& vertex = id_vertex.second;
      auto local_scope = vertex_ids_with_local_scope.find(vid);

      if (local_scope != vertex_ids_with_local_scope.end())
        grid.vertices.Insert(vid, Vector3(vertex.x, vertex.y, layer_z_level));

      ++vid;
    } // for vertex
  }   // for layer

  grid.SetGlobalVertexCount(vid);
}

void
VolumeMesherExtruder::Execute()
{
  Chi::log.Log() << Chi::program_timer.GetTimeString()
                 << " VolumeMesherExtruder executed. Memory in use = " << Chi::GetMemoryUsageInMB()
                 << " MB" << std::endl;

  // Loop over all regions
  Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Processing Region" << std::endl;

  // Create new continuum
  auto grid = MeshContinuum::New();
  auto temp_grid = MeshContinuum::New();

  SetContinuum(grid);
  SetGridAttributes(DIMENSION_3 | EXTRUDED);

  // Setup layers
  // populates vertex-layers
  Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Setting up layers" << std::endl;
  SetupLayers();

  // Process templates
  if (template_type_ == TemplateType::UNPARTITIONED_MESH)
  {
    Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Processing unpartitioned mesh" << std::endl;

    // Get node_z_incr
    node_z_index_incr_ = template_unpartitioned_mesh_->GetVertices().size();

    // Create baseline polygons in template continuum
    Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Creating template cells" << std::endl;
    CreatePolygonCells(*template_unpartitioned_mesh_, temp_grid);

    grid->GetBoundaryIDMap() = template_unpartitioned_mesh_->GetMeshOptions().boundary_id_map;
    temp_grid->GetBoundaryIDMap() = template_unpartitioned_mesh_->GetMeshOptions().boundary_id_map;
  }

  Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Creating local nodes" << std::endl;
  CreateLocalNodes(*temp_grid, *grid);

  Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Done creating local nodes" << std::endl;

  // Insert top and bottom boundary id map
  auto& grid_bndry_id_map = grid->GetBoundaryIDMap();
  zmax_bndry_id = grid->MakeBoundaryID("ZMAX");
  grid_bndry_id_map[zmax_bndry_id] = "ZMAX";
  zmin_bndry_id = grid->MakeBoundaryID("ZMIN");
  grid_bndry_id_map[zmin_bndry_id] = "ZMIN";

  // Create extruded item_id
  Chi::log.Log() << "VolumeMesherExtruder: Extruding cells" << std::endl;
  opensn::mpi.Barrier();
  ExtrudeCells(*temp_grid, *grid);

  size_t total_local_cells = grid->local_cells.size();
  size_t total_global_cells = 0;

  MPI_Allreduce(
    &total_local_cells, &total_global_cells, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi.comm);

  Chi::log.Log() << "VolumeMesherExtruder: Cells extruded = " << total_global_cells << std::endl;

  // Checking partitioning parameters
  if (options.partition_type != KBA_STYLE_XYZ)
  {
    Chi::log.LogAllError() << "Any partitioning scheme other than KBA_STYLE_XYZ is currently not"
                              " supported by VolumeMesherExtruder. No worries. There are plans"
                              " to develop this support.";
    Chi::Exit(EXIT_FAILURE);
  }
  if (!options.mesh_global)
  {
    int p_tot = options.partition_x * options.partition_y * options.partition_z;

    if (opensn::mpi.process_count != p_tot)
    {
      Chi::log.LogAllError() << "ERROR: Number of processors available ("
                             << opensn::mpi.process_count
                             << ") does not match amount of processors "
                             << "required by surface mesher partitioning parameters (" << p_tot
                             << ").";
      Chi::Exit(EXIT_FAILURE);
    }
  } // if mesh-global

  Chi::log.LogAllVerbose1() << "Building local cell indices";

  // Print info
  Chi::log.LogAllVerbose1() << "### LOCATION[" << opensn::mpi.location_id
                            << "] amount of local cells=" << grid->local_cells.size();

  Chi::log.Log() << "VolumeMesherExtruder: Number of cells in region = " << total_global_cells
                 << std::endl;

  opensn::mpi.Barrier();
}

void
VolumeMesherExtruder::ExtrudeCells(MeshContinuum& template_grid, MeshContinuum& grid)
{
  const Vector3 khat(0.0, 0.0, 1.0);
  // Start extrusion
  size_t num_global_cells = 0;
  for (size_t iz = 0; iz < (vertex_layers_.size() - 1); iz++)
  {
    for (const auto& template_cell : template_grid.local_cells)
    {
      // Check template cell type
      if (template_cell.Type() != CellType::POLYGON)
        throw std::logic_error("Extruder::ExtrudeCells: "
                               "Template cell error. Not of base type POLYGON");

      // Check cell not inverted
      {
        const auto& v0 = template_cell.centroid_;
        const auto& v1 = template_grid.vertices[template_cell.vertex_ids_[0]];
        const auto& v2 = template_grid.vertices[template_cell.vertex_ids_[1]];

        auto v01 = v1 - v0;
        auto v02 = v2 - v0;

        if (v01.Cross(v02).Dot(khat) < 0.0)
          throw std::logic_error("Extruder attempting to extrude a template"
                                 " cell with a normal pointing downward. This"
                                 " causes erratic behavior and needs to be"
                                 " corrected.");
      }

      auto projected_centroid = ProjectCentroidToLevel(template_cell.centroid_, iz);
      int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);

      bool has_local_scope = HasLocalScope(template_cell, template_grid, iz);

      if (has_local_scope)
      {
        auto cell = MakeExtrudedCell(
          template_cell, grid, iz, num_global_cells, pid, template_grid.local_cells.size());

        cell->material_id_ = template_cell.material_id_;

        grid.cells.push_back(std::move(cell));
      }
      ++num_global_cells;

    } // for template cell

  } // for iz
}

void
VolumeMesherExtruder::SetupLayers(int default_layer_count)
{
  // Create default layers if no input layers are provided
  if (input_layers_.empty())
  {
    Chi::log.Log0Warning() << "VolumeMesherExtruder: No extrusion layers have been specified. "
                           << "A default single layer will be used with height 1.0 and a single "
                           << "subdivision.";
    double dz = 1.0 / default_layer_count;
    for (int k = 0; k <= default_layer_count; k++)
    {
      vertex_layers_.push_back(k * dz);
    }
  }
  else
  {
    double last_z = 0.0;
    vertex_layers_.push_back(last_z);

    for (const auto& input_layer : input_layers_)
    {
      double dz = input_layer.height / input_layer.sub_divisions;

      for (int k = 0; k < input_layer.sub_divisions; k++)
      {
        last_z += dz;
        vertex_layers_.push_back(last_z);
      }
    }
  }

  Chi::log.Log() << "VolumeMesherExtruder: Total number of cell layers is "
                 << vertex_layers_.size() - 1;
}

Vector3
VolumeMesherExtruder::ProjectCentroidToLevel(const Vector3& centroid, const size_t level)
{
  double z_location = 0.5 * (vertex_layers_[level] + vertex_layers_[level + 1]);

  auto centroid_projected = centroid;
  centroid_projected.z = z_location;

  return centroid_projected;
}

int
VolumeMesherExtruder::GetCellKBAPartitionIDFromCentroid(Vector3& centroid)
{
  int px = options.partition_x;
  int py = options.partition_y;

  Cell n_gcell(CellType::GHOST, CellType::GHOST);
  n_gcell.centroid_ = centroid;

  auto xyz_partition_indices = GetCellXYZPartitionID(&n_gcell);

  int nxi = std::get<0>(xyz_partition_indices);
  int nyi = std::get<1>(xyz_partition_indices);
  int nzi = std::get<2>(xyz_partition_indices);

  return nzi * px * py + nyi * px + nxi;
}

bool
VolumeMesherExtruder::HasLocalScope(const Cell& template_cell,
                                    const MeshContinuum& template_continuum,
                                    size_t z_level)
{
  // Check if the template cell is in the current partition
  {
    auto& centroid = template_cell.centroid_;
    auto projected_centroid = ProjectCentroidToLevel(centroid, z_level);
    int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
    if (pid == opensn::mpi.location_id) return true;
  }

  const size_t last_z_level = vertex_layers_.size() - 1;
  const size_t z_level_below = z_level - 1;
  const size_t z_level_above = z_level + 1;

  // Build z-levels to search
  std::vector<size_t> z_levels_to_search;
  z_levels_to_search.reserve(3);

  z_levels_to_search.push_back(z_level);
  if (z_level != 0) z_levels_to_search.push_back(z_level_below);
  if (z_level != (last_z_level - 1)) z_levels_to_search.push_back(z_level_above);

  // Search template cell's lateral neighbors
  const auto& vertex_subs = template_unpartitioned_mesh_->GetVertextCellSubscriptions();
  for (uint64_t vid : template_cell.vertex_ids_)
    for (uint64_t cid : vertex_subs[vid])
    {
      if (cid == template_cell.local_id_) continue;

      auto& candidate_cell = template_continuum.local_cells[cid];
      auto& cc_centroid = candidate_cell.centroid_;

      for (size_t z : z_levels_to_search)
      {
        auto projected_centroid = ProjectCentroidToLevel(cc_centroid, z);
        int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
        if (pid == opensn::mpi.location_id) return true;
      }
    } // for cid

  // Search template cell's longitudinal neighbors
  for (size_t z : z_levels_to_search)
  {
    if (z == z_level) continue;

    auto projected_centroid = ProjectCentroidToLevel(template_cell.centroid_, z);
    int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);
    if (pid == opensn::mpi.location_id) return true;
  } // for z

  return false;
}

std::unique_ptr<Cell>
VolumeMesherExtruder::MakeExtrudedCell(const Cell& template_cell,
                                       const MeshContinuum& grid,
                                       size_t z_level,
                                       uint64_t cell_global_id,
                                       int partition_id,
                                       size_t num_template_cells)
{
  const size_t tc_num_verts = template_cell.vertex_ids_.size();

  // Determine cell sub-type
  CellType extruded_subtype;
  switch (template_cell.SubType())
  {
    case CellType::TRIANGLE:
      extruded_subtype = CellType::WEDGE;
      break;
    case CellType::QUADRILATERAL:
      extruded_subtype = CellType::HEXAHEDRON;
      break;
    default:
      extruded_subtype = CellType::POLYHEDRON;
  }

  // Create polyhedron
  auto cell = std::make_unique<Cell>(CellType::POLYHEDRON, extruded_subtype);
  cell->global_id_ = cell_global_id;
  // cell->local_id set when added to mesh
  cell->partition_id_ = partition_id;
  cell->centroid_ = ProjectCentroidToLevel(template_cell.centroid_, z_level);

  // Populate cell v-indices
  cell->vertex_ids_.reserve(2 * tc_num_verts);
  for (auto tc_vid : template_cell.vertex_ids_)
    cell->vertex_ids_.push_back(tc_vid + z_level * node_z_index_incr_);

  for (auto tc_vid : template_cell.vertex_ids_)
    cell->vertex_ids_.push_back(tc_vid + (z_level + 1) * node_z_index_incr_);

  // Create side faces
  for (auto& face : template_cell.faces_)
  {
    CellFace newFace;

    newFace.vertex_ids_.resize(4, -1);
    newFace.vertex_ids_[0] = face.vertex_ids_[0] + z_level * node_z_index_incr_;
    newFace.vertex_ids_[1] = face.vertex_ids_[1] + z_level * node_z_index_incr_;
    newFace.vertex_ids_[2] = face.vertex_ids_[1] + (z_level + 1) * node_z_index_incr_;
    newFace.vertex_ids_[3] = face.vertex_ids_[0] + (z_level + 1) * node_z_index_incr_;

    // Compute centroid
    const auto& v0 = grid.vertices[newFace.vertex_ids_[0]];
    const auto& v1 = grid.vertices[newFace.vertex_ids_[1]];
    const auto& v2 = grid.vertices[newFace.vertex_ids_[2]];
    const auto& v3 = grid.vertices[newFace.vertex_ids_[3]];

    Vertex vfc = (v0 + v1 + v2 + v3) / 4.0;
    newFace.centroid_ = vfc;

    // Compute normal
    Vector3 va = v0 - vfc;
    Vector3 vb = v1 - vfc;

    Vector3 vn = va.Cross(vb);

    newFace.normal_ = (vn / vn.Norm());

    // Set neighbor
    // The side connections have the same connections as the
    // template cell + the z_level specifiers of the layer.
    if (face.has_neighbor_)
    {
      newFace.neighbor_id_ = face.neighbor_id_ + z_level * num_template_cells;
      newFace.has_neighbor_ = true;
    }
    else
      newFace.neighbor_id_ = face.neighbor_id_;

    cell->faces_.push_back(newFace);
  } // for side faces

  // Create top and bottom faces
  // Bottom face
  {
    CellFace newFace;
    // Vertices
    auto vfc = Vertex(0.0, 0.0, 0.0);
    newFace.vertex_ids_.reserve(tc_num_verts);
    for (int tv = (static_cast<int>(tc_num_verts) - 1); tv >= 0; tv--)
    {
      newFace.vertex_ids_.push_back(template_cell.vertex_ids_[tv] + z_level * node_z_index_incr_);
      const auto& v = grid.vertices[newFace.vertex_ids_.back()];
      vfc = vfc + v;
    }

    // Compute centroid
    vfc = vfc / static_cast<double>(tc_num_verts);
    newFace.centroid_ = vfc;

    // Compute normal
    auto va = grid.vertices[newFace.vertex_ids_[0]] - vfc;
    auto vb = grid.vertices[newFace.vertex_ids_[1]] - vfc;

    auto vn = va.Cross(vb);
    newFace.normal_ = vn / vn.Norm();

    // Set neighbor
    if (z_level == 0)
    {
      newFace.neighbor_id_ = zmin_bndry_id;
      newFace.has_neighbor_ = false;
    }
    else
    {
      newFace.neighbor_id_ = template_cell.local_id_ + (z_level - 1) * num_template_cells;
      newFace.has_neighbor_ = true;
    }

    cell->faces_.push_back(newFace);
  }

  // Top face
  {
    CellFace newFace;
    // Vertices
    auto vfc = Vertex(0.0, 0.0, 0.0);
    newFace.vertex_ids_.reserve(tc_num_verts);
    for (auto tc_vid : template_cell.vertex_ids_)
    {
      newFace.vertex_ids_.push_back(tc_vid + (z_level + 1) * node_z_index_incr_);
      const auto& v = grid.vertices[newFace.vertex_ids_.back()];
      vfc = vfc + v;
    }

    // Compute centroid
    vfc = vfc / static_cast<double>(tc_num_verts);
    newFace.centroid_ = vfc;

    // Compute normal
    auto va = grid.vertices[newFace.vertex_ids_[0]] - vfc;
    auto vb = grid.vertices[newFace.vertex_ids_[1]] - vfc;

    auto vn = va.Cross(vb);
    newFace.normal_ = vn / vn.Norm();

    // Set neighbor
    if (z_level == (vertex_layers_.size() - 2))
    {
      newFace.neighbor_id_ = zmax_bndry_id;
      newFace.has_neighbor_ = false;
    }
    else
    {
      newFace.neighbor_id_ = template_cell.local_id_ + (z_level + 1) * num_template_cells;
      newFace.has_neighbor_ = true;
    }

    cell->faces_.push_back(newFace);
  }

  return cell;
}

} // namespace opensn
