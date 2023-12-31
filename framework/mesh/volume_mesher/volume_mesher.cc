#include "framework/mesh/volume_mesher/volume_mesher.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/volume_mesher/extruder/volmesher_extruder.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/console/console.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi.h"

namespace opensn
{

VolumeMesher::VolumeMesher(VolumeMesherType type) : type_(type)
{
}

void
VolumeMesher::SetContinuum(MeshContinuumPtr& grid)
{
  grid_ptr_ = grid;
}

MeshContinuumPtr&
VolumeMesher::GetContinuum()
{
  return grid_ptr_;
}

void
VolumeMesher::SetGridAttributes(MeshAttributes new_attribs, std::array<size_t, 3> ortho_Nis)
{
  grid_ptr_->SetAttributes(new_attribs, ortho_Nis);
}

VolumeMesherType
VolumeMesher::Type() const
{
  return type_;
}

std::pair<int, int>
VolumeMesher::GetCellXYPartitionID(Cell* cell)
{
  std::pair<int, int> ij_id(0, 0);

  if (opensn::mpi.process_count == 1) { return ij_id; }

  // Get the current handler
  auto& mesh_handler = GetCurrentHandler();
  auto& vol_mesher = mesh_handler.GetVolumeMesher();

  // Sanity check on partitioning
  size_t num_x_subsets = vol_mesher.options.xcuts.size() + 1;
  size_t num_y_subsets = vol_mesher.options.ycuts.size() + 1;

  size_t x_remainder = num_x_subsets % vol_mesher.options.partition_x;
  size_t y_remainder = num_y_subsets % vol_mesher.options.partition_y;

  if (x_remainder != 0)
  {
    log.LogAllError() << "When specifying x-partitioning, the number of grp_subsets in x "
                         "needs to be divisible by the number of partitions in x.";
    Exit(EXIT_FAILURE);
  }

  if (y_remainder != 0)
  {
    log.LogAllError() << "When specifying y-partitioning, the number of grp_subsets in y "
                         "needs to be divisible by the number of partitions in y.";
    Exit(EXIT_FAILURE);
  }

  size_t subsets_per_partitionx = num_x_subsets / vol_mesher.options.partition_x;
  size_t subsets_per_partitiony = num_y_subsets / vol_mesher.options.partition_y;

  // Determine x-partition
  int x = -1;
  int xcount = -1;
  for (size_t i = subsets_per_partitionx - 1; i < vol_mesher.options.xcuts.size();
       i += subsets_per_partitionx)
  {
    xcount++;
    if (cell->centroid_.x <= vol_mesher.options.xcuts[i])
    {
      x = xcount;
      break;
    }
  }
  if (x < 0) { x = vol_mesher.options.partition_x - 1; }

  // Determine y-partition
  int y = -1;
  int ycount = -1;
  for (size_t i = subsets_per_partitiony - 1; i < vol_mesher.options.ycuts.size();
       i += subsets_per_partitiony)
  {
    ycount++;
    if (cell->centroid_.y <= vol_mesher.options.ycuts[i])
    {
      y = ycount;
      break;
    }
  }
  if (y < 0) { y = vol_mesher.options.partition_y - 1; }

  // Set partitioning
  ij_id.first = x;
  ij_id.second = y;

  return ij_id;
}

std::tuple<int, int, int>
VolumeMesher::GetCellXYZPartitionID(Cell* cell)
{
  std::tuple<int, int, int> ijk_id(0, 0, 0);
  bool found_partition = false;

  if (opensn::mpi.process_count == 1) { return ijk_id; }

  // Get ij indices
  std::pair<int, int> ij_id = GetCellXYPartitionID(cell);

  // Get the current handler
  auto& mesh_handler = GetCurrentHandler();
  auto& vol_mesher = mesh_handler.GetVolumeMesher();

  if (vol_mesher.options.partition_z == 1)
  {
    found_partition = true;
    std::get<0>(ijk_id) = ij_id.first;
    std::get<1>(ijk_id) = ij_id.second;
    std::get<2>(ijk_id) = 0;
  }
  else if (vol_mesher.Type() == VolumeMesherType::EXTRUDER)
  {
    auto extruder = dynamic_cast<VolumeMesherExtruder&>(vol_mesher);
    const auto& vertex_layers = extruder.GetVertexLayers();
    // Create virtual cuts
    if (vol_mesher.options.zcuts.empty())
    {
      size_t num_sub_layers = vertex_layers.size() - 1;

      if ((num_sub_layers % vol_mesher.options.partition_z) != 0)
      {
        log.LogAllError() << "Number of sub-layers in extruded mesh is not divisible "
                          << "by the requested number of z-partitions.";
        Exit(EXIT_FAILURE);
      }

      int delta_zk = num_sub_layers / vol_mesher.options.partition_z;
      for (int k = 0; k < (vol_mesher.options.partition_z); k++)
      {
        int layer_index = k * delta_zk + delta_zk;
        if (layer_index > (vertex_layers.size() - 1))
        {
          layer_index = (int)vertex_layers.size() - 1;
          vol_mesher.options.zcuts.push_back(vertex_layers[layer_index]);
        }
        else
        {
          vol_mesher.options.zcuts.push_back(vertex_layers[layer_index]);

          if (log.GetVerbosity() == Logger::LOG_LVL::LOG_0VERBOSE_2)
          {
            printf("Z-Cut %lu, %g\n", vol_mesher.options.zcuts.size(), vertex_layers[layer_index]);
          }
        }
      }
    }

    // Scan cuts for location
    double zmin = -1.0e-16;
    for (int k = 0; k < (vol_mesher.options.zcuts.size()); k++)
    {
      double zmax = vol_mesher.options.zcuts[k];

      double z = cell->centroid_.z;

      if (log.GetVerbosity() == Logger::LOG_0VERBOSE_2)
      {
        printf("zmax = %g, zmin = %g, cell_z = %g\n", zmax, zmin, z);
      }

      if ((z > zmin) && (z < zmax))
      {
        std::get<0>(ijk_id) = ij_id.first;
        std::get<1>(ijk_id) = ij_id.second;
        std::get<2>(ijk_id) = k;

        found_partition = true;
        break;
      }
      zmin = zmax;
    }
  } // if typeid
  else if (vol_mesher.Type() == VolumeMesherType::UNPARTITIONED)
  {
    if (vol_mesher.options.zcuts.empty())
    {
      throw std::invalid_argument("Cell z-partitioning cannot be determined "
                                  "because no z-cuts are supplied to volume "
                                  "mesher.");
    }

    // Scan cuts for location
    std::vector<double> temp_zcuts = vol_mesher.options.zcuts;
    double zmin = -1.0e16;
    double zmax = 1.0e16;
    temp_zcuts.push_back(zmax);
    for (int k = 0; k < (temp_zcuts.size()); k++)
    {
      zmax = temp_zcuts[k];

      double z = cell->centroid_.z;

      if (log.GetVerbosity() == Logger::LOG_0VERBOSE_2)
      {
        printf("zmax = %g, zmin = %g, cell_z = %g\n", zmax, zmin, z);
      }

      if ((z > zmin) && (z < zmax))
      {
        std::get<0>(ijk_id) = ij_id.first;
        std::get<1>(ijk_id) = ij_id.second;
        std::get<2>(ijk_id) = k;

        found_partition = true;
        break;
      }
      zmin = zmax;
    } // for k
  }   // if typeid

  // Report unallocated
  // item_id
  if (!found_partition)
  {
    log.LogAllError() << "A cell was encountered for which "
                         "no zpartition id was found";
    Exit(EXIT_FAILURE);
  }

  return ijk_id;
}

void
VolumeMesher::CreatePolygonCells(const UnpartitionedMesh& umesh, MeshContinuumPtr& grid)
{
  // Copy nodes
  {
    uint64_t id = 0;
    for (const auto& vertex : umesh.GetVertices())
      grid->vertices.Insert(id++, vertex);
  }

  size_t num_cells = 0;
  for (auto& raw_cell : umesh.GetRawCells())
  {
    // Check valid template cell
    if (raw_cell->type != CellType::POLYGON)
    {
      log.LogAllError() << "VolumeMesher::CreatePolygonCells "
                           "called with a cell not being of primary type"
                           " CellType::POLYGON.";
      Exit(EXIT_FAILURE);
    }

    // Make cell
    auto cell = std::make_unique<Cell>(CellType::POLYGON, raw_cell->sub_type);

    cell->global_id_ = num_cells;
    cell->local_id_ = num_cells;
    cell->partition_id_ = opensn::mpi.location_id;

    cell->centroid_ = raw_cell->centroid;
    cell->material_id_ = raw_cell->material_id;
    cell->vertex_ids_ = raw_cell->vertex_ids;

    // Copy faces + compute face centroid and normal
    const Vector3 khat(0.0, 0.0, 1.0);
    for (auto& raw_face : raw_cell->faces)
    {
      CellFace new_face;

      new_face.vertex_ids_ = raw_face.vertex_ids;

      const auto& v0 = grid->vertices[new_face.vertex_ids_[0]];
      const auto& v1 = grid->vertices[new_face.vertex_ids_[1]];
      new_face.centroid_ = v0 * 0.5 + v1 * 0.5;

      Vector3 va = v1 - v0;
      Vector3 vn = va.Cross(khat);
      vn = vn / vn.Norm();
      new_face.normal_ = vn;

      new_face.has_neighbor_ = raw_face.has_neighbor;
      new_face.neighbor_id_ = raw_face.neighbor;

      cell->faces_.push_back(new_face);
    }

    // Push to grid
    grid->cells.push_back(std::move(cell));
    ++num_cells;
  } // for raw_cell
}

void
VolumeMesher::SetMatIDFromLogical(const LogicalVolume& log_vol, bool sense, int mat_id)
{
  log.Log0Verbose1() << program_timer.GetTimeString()
                     << " Setting material id from logical volume.";
  // Get current mesh handler
  auto& handler = GetCurrentHandler();

  // Get back mesh
  MeshContinuumPtr vol_cont = handler.GetGrid();

  int num_cells_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    if (log_vol.Inside(cell.centroid_) && sense)
    {
      cell.material_id_ = mat_id;
      ++num_cells_modified;
    }
  }

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = vol_cont->cells[ghost_id];
    if (log_vol.Inside(cell.centroid_) && sense) cell.material_id_ = mat_id;
  }

  int global_num_cells_modified;
  MPI_Allreduce(&num_cells_modified, &global_num_cells_modified, 1, MPI_INT, MPI_SUM, mpi.comm);

  log.Log0Verbose1() << program_timer.GetTimeString()
                     << " Done setting material id from logical volume. "
                     << "Number of cells modified = " << global_num_cells_modified << ".";
}

void
VolumeMesher::SetBndryIDFromLogical(const LogicalVolume& log_vol,
                                    bool sense,
                                    const std::string& bndry_name)
{
  log.Log() << program_timer.GetTimeString() << " Setting boundary id from logical volume.";
  // Get current mesh handler
  auto& handler = GetCurrentHandler();

  // Get back mesh
  MeshContinuumPtr vol_cont = handler.GetGrid();

  // Check if name already has id
  auto& grid_bndry_id_map = vol_cont->GetBoundaryIDMap();
  uint64_t bndry_id = vol_cont->MakeBoundaryID(bndry_name);

  // Loop over cells
  int num_faces_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    for (auto& face : cell.faces_)
    {
      if (face.has_neighbor_) continue;
      if (log_vol.Inside(face.centroid_) && sense)
      {
        face.neighbor_id_ = bndry_id;
        ++num_faces_modified;
      }
    }
  }

  int global_num_faces_modified;
  MPI_Allreduce(&num_faces_modified, &global_num_faces_modified, 1, MPI_INT, MPI_SUM, mpi.comm);

  if (global_num_faces_modified > 0 and grid_bndry_id_map.count(bndry_id) == 0)
    grid_bndry_id_map[bndry_id] = bndry_name;

  log.Log() << program_timer.GetTimeString() << " Done setting boundary id from logical volume. "
            << "Number of faces modified = " << global_num_faces_modified << ".";
}

void
VolumeMesher::SetMatIDToAll(int mat_id)
{
  log.Log() << program_timer.GetTimeString() << " Setting material id " << mat_id
            << " to all cells.";

  // Get current mesh handler
  auto& handler = GetCurrentHandler();

  // Get back mesh
  auto vol_cont = handler.GetGrid();

  for (auto& cell : vol_cont->local_cells)
    cell.material_id_ = mat_id;

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    vol_cont->cells[ghost_id].material_id_ = mat_id;

  opensn::mpi.Barrier();
  log.Log() << program_timer.GetTimeString() << " Done setting material id " << mat_id
            << " to all cells";
}

void
VolumeMesher::SetupOrthogonalBoundaries()
{
  log.Log() << program_timer.GetTimeString() << " Setting orthogonal boundaries.";

  // Get current mesh handler
  auto& handler = GetCurrentHandler();

  // Get back mesh
  auto vol_cont = handler.GetGrid();

  const Vector3 ihat(1.0, 0.0, 0.0);
  const Vector3 jhat(0.0, 1.0, 0.0);
  const Vector3 khat(0.0, 0.0, 1.0);

  for (auto& cell : vol_cont->local_cells)
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        Vector3& n = face.normal_;

        std::string boundary_name;
        if (n.Dot(ihat) > 0.999) boundary_name = "XMAX";
        else if (n.Dot(ihat) < -0.999)
          boundary_name = "XMIN";
        else if (n.Dot(jhat) > 0.999)
          boundary_name = "YMAX";
        else if (n.Dot(jhat) < -0.999)
          boundary_name = "YMIN";
        else if (n.Dot(khat) > 0.999)
          boundary_name = "ZMAX";
        else if (n.Dot(khat) < -0.999)
          boundary_name = "ZMIN";

        uint64_t bndry_id = vol_cont->MakeBoundaryID(boundary_name);

        face.neighbor_id_ = bndry_id;

        vol_cont->GetBoundaryIDMap()[bndry_id] = boundary_name;
      } // if bndry

  opensn::mpi.Barrier();
  log.Log() << program_timer.GetTimeString() << " Done setting orthogonal boundaries.";
}

void
VolumeMesher::Execute()
{
  log.Log() << "Empty volume mesher, nothing executed.";
  log.Log() << std::endl;
}

} // namespace opensn
