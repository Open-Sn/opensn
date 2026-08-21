// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_generator/distributed_mesh_generator.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/data_types/byte_array.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"

namespace opensn
{

DistributedMeshGenerator::DistributedMeshGenerator(const InputParameters& params)
  : MeshGenerator(params), num_partitions_(mpi_comm.size())
{
}

std::shared_ptr<Mesh>
DistributedMeshGenerator::Execute()
{
  const auto rank = mpi_comm.rank();
  const auto num_partitions = mpi_comm.size();
  DistributedMeshData mesh_info;

  log.Log() << program_timer.GetTimeString() << " Distributing mesh with " << num_partitions
            << " parts";

  if (rank == 0)
  {
    std::shared_ptr<UnpartitionedMesh> current_umesh = nullptr;
    for (const auto& mesh_generator_ptr : inputs_)
      current_umesh = mesh_generator_ptr->GenerateUnpartitionedMesh(current_umesh);
    current_umesh = GenerateUnpartitionedMesh(current_umesh);

    const auto cell_pids = PartitionMesh(*current_umesh, num_partitions);
    auto serial_data = DistributeSerializedMeshData(cell_pids, *current_umesh, num_partitions);
    mesh_info = DeserializeMeshData(serial_data);
  }
  else
  {
    std::vector<std::byte> data;
    mpi_comm.recv<std::byte>(0, rank, data);
    ByteArray serial_data(data);
    mesh_info = DeserializeMeshData(serial_data);
  }

  auto grid_ptr = SetupLocalMesh(mesh_info);

  mpi_comm.barrier();

  log.Log() << program_timer.GetTimeString() << " Mesh successfully distributed";

  return grid_ptr;
}

OpenSnRegisterObjectInNamespace(mesh, DistributedMeshGenerator);

InputParameters
DistributedMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription(
    "Generates and partitions the mesh on location 0. The partitioned mesh is "
    "broadcast to all other locations.");

  return params;
}

std::shared_ptr<DistributedMeshGenerator>
DistributedMeshGenerator::Create(const ParameterBlock& params)
{
  const auto& factory = ObjectFactory::GetInstance();
  return factory.Create<DistributedMeshGenerator>("mesh::DistributedMeshGenerator", params);
}

ByteArray
DistributedMeshGenerator::DistributeSerializedMeshData(const std::vector<int>& cell_pids,
                                                       const UnpartitionedMesh& umesh,
                                                       const int num_partitions)
{
  const auto& vertex_subs = umesh.GetVertextCellSubscriptions();
  const auto& raw_cells = umesh.GetCells();
  const auto& raw_vertices = umesh.GetVertices();
  ByteArray loc0_data;

  for (int pid = 0; pid < num_partitions; ++pid)
  {
    ByteArray serial_data;

    std::vector<uint64_t> local_cells_needed;
    std::set<uint64_t> cells_needed;
    std::set<uint64_t> vertices_needed;
    local_cells_needed.reserve(cell_pids.size() / num_partitions);

    for (uint64_t cell_global_id = 0; cell_global_id < cell_pids.size(); ++cell_global_id)
    {
      if (cell_pids[cell_global_id] == pid)
      {
        const auto& raw_cell = raw_cells[cell_global_id];
        local_cells_needed.push_back(cell_global_id);
        cells_needed.emplace(cell_global_id);

        auto raw_cell_vertex_ids = umesh.GetCellConnectivity(cell_global_id);
        for (const auto vid : raw_cell_vertex_ids)
        {
          vertices_needed.emplace(vid);

          // Process ghost cells
          for (const auto ghost_gid : vertex_subs[vid])
          {
            if (ghost_gid != cell_global_id && cells_needed.find(ghost_gid) == cells_needed.end())
            {
              cells_needed.emplace(ghost_gid);
              auto ghost_raw_cell_vertex_ids = umesh.GetCellConnectivity(ghost_gid);

              // Insert ghost vertex IDs
              for (const auto gvid : ghost_raw_cell_vertex_ids)
                vertices_needed.emplace(gvid);
            }
          }
        }
      }
    }

    // Basic mesh data
    serial_data.Write<unsigned int>(umesh.GetDimension());
    serial_data.Write(static_cast<int>(umesh.GetCoordinateSystem()));
    serial_data.Write(static_cast<int>(umesh.GetType()));
    serial_data.Write(umesh.IsExtruded());
    const auto& [Nx, Ny, Nz] = umesh.GetOrthoAttributes();
    serial_data.Write(Nx);
    serial_data.Write(Ny);
    serial_data.Write(Nz);
    serial_data.Write(raw_vertices.size());

    // Boundaries
    const auto& bndry_map = umesh.GetBoundaryIDMap();
    serial_data.Write(bndry_map.size());
    for (const auto& [bid, bname] : bndry_map)
    {
      serial_data.Write(bid);
      const size_t num_chars = bname.size();
      serial_data.Write(num_chars);
      for (size_t i = 0; i < num_chars; ++i)
        serial_data.Write(bname[i]);
    }

    // Number of cells and vertices
    serial_data.Write(cells_needed.size());
    serial_data.Write(vertices_needed.size());

    // Cell data
    for (const auto& cell_global_id : cells_needed)
    {
      const auto& cell = raw_cells[cell_global_id];
      serial_data.Write(static_cast<int>(cell_pids[cell_global_id]));
      serial_data.Write(cell_global_id);
      serial_data.Write(cell.GetType());
      serial_data.Write(cell.GetSubType());
      serial_data.Write(cell.centroid.x);
      serial_data.Write(cell.centroid.y);
      serial_data.Write(cell.centroid.z);
      serial_data.Write(cell.block_id);

      auto cell_vertex_ids = umesh.GetCellConnectivity(cell_global_id);
      serial_data.Write(cell_vertex_ids.size());
      for (const auto vid : cell_vertex_ids)
        serial_data.Write(vid);

      const auto num_faces = umesh.GetCellFaceCount(cell_global_id);
      serial_data.Write(num_faces);
      for (std::size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = umesh.GetCellFace(cell_global_id, f);
        const auto& face_vids = umesh.GetCellFaceConnectivity()[cell_global_id][f];
        serial_data.Write(face_vids.size());
        for (const auto vid : face_vids)
          serial_data.Write(vid);
        serial_data.Write(face.has_neighbor);
        serial_data.Write(face.neighbor_id);
      }
    }

    // Vertex data
    for (const auto vid : vertices_needed)
    {
      serial_data.Write(vid);
      serial_data.Write(raw_vertices[vid].x);
      serial_data.Write(raw_vertices[vid].y);
      serial_data.Write(raw_vertices[vid].z);
    }

    if (pid == 0)
      loc0_data = serial_data;
    else
      mpi_comm.send<std::byte>(
        pid, pid, serial_data.Data().data(), static_cast<int>(serial_data.Size()));
  }
  return loc0_data;
}

DistributedMeshGenerator::DistributedMeshData
DistributedMeshGenerator::DeserializeMeshData(ByteArray& serial_data)
{
  DistributedMeshData info_block;

  // Basic mesh data
  info_block.dimension = serial_data.Read<unsigned int>();
  info_block.coord_sys = static_cast<CoordinateSystemType>(serial_data.Read<int>());
  info_block.mesh_type = static_cast<MeshType>(serial_data.Read<int>());
  info_block.extruded = serial_data.Read<bool>();
  info_block.ortho_attributes.Nx = serial_data.Read<size_t>();
  info_block.ortho_attributes.Ny = serial_data.Read<size_t>();
  info_block.ortho_attributes.Nz = serial_data.Read<size_t>();
  info_block.num_global_vertices = serial_data.Read<size_t>();

  // Boundaries
  auto num_boundaries = serial_data.Read<size_t>();
  for (size_t b = 0; b < num_boundaries; ++b)
  {
    const auto bid = serial_data.Read<uint64_t>();
    const auto num_chars = serial_data.Read<size_t>();
    std::string bname(num_chars, ' ');
    for (size_t i = 0; i < num_chars; ++i)
      bname[i] = serial_data.Read<char>();
    info_block.boundary_id_map.insert(std::make_pair(bid, bname));
  }

  // Number of cells and vertices
  const auto num_cells = serial_data.Read<size_t>();
  const auto num_vertices = serial_data.Read<size_t>();

  // Cell data
  for (size_t i = 0; i < num_cells; ++i)
  {
    const auto cell_pid = serial_data.Read<int>();
    const auto cell_gid = serial_data.Read<uint64_t>();
    const auto type = serial_data.Read<CellType>();
    const auto sub_type = serial_data.Read<CellType>();

    Cell cell(type, sub_type);

    cell.centroid.x = serial_data.Read<double>();
    cell.centroid.y = serial_data.Read<double>();
    cell.centroid.z = serial_data.Read<double>();
    cell.block_id = serial_data.Read<unsigned int>();

    const auto num_vids = serial_data.Read<size_t>();
    std::vector<std::uint64_t> cell_vertex_ids;
    cell_vertex_ids.reserve(num_vids);
    for (size_t v = 0; v < num_vids; ++v)
      cell_vertex_ids.push_back(serial_data.Read<uint64_t>());
    info_block.cell_connect.emplace(cell_gid, cell_vertex_ids);

    const auto num_faces = serial_data.Read<size_t>();
    std::vector<std::vector<std::uint64_t>> cell_faces_vids;
    cell_faces_vids.reserve(num_faces);
    std::vector<CellFace> cell_faces;
    cell_faces.reserve(num_faces);
    for (size_t f = 0; f < num_faces; ++f)
    {
      CellFace face;
      auto num_face_vids = serial_data.Read<size_t>();
      std::vector<std::uint64_t> face_vids;
      face_vids.reserve(num_face_vids);
      for (size_t v = 0; v < num_face_vids; ++v)
        face_vids.push_back(serial_data.Read<uint64_t>());
      cell_faces_vids.push_back(std::move(face_vids));

      face.has_neighbor = serial_data.Read<bool>();
      face.neighbor_id = serial_data.Read<uint64_t>();

      cell_faces.emplace_back(face);
    }
    info_block.cell_face_connect.emplace(cell_gid, cell_faces_vids);
    info_block.cell_faces.emplace(cell_gid, std::move(cell_faces));
    info_block.cells.insert(std::make_pair(std::make_pair(cell_pid, cell_gid), cell));
  }

  // Vertex data
  for (size_t i = 0; i < num_vertices; ++i)
  {
    auto vid = serial_data.Read<uint64_t>();
    Vector3 vertex;
    vertex.x = serial_data.Read<double>();
    vertex.y = serial_data.Read<double>();
    vertex.z = serial_data.Read<double>();
    info_block.vertices.insert(std::make_pair(vid, vertex));
  }
  return info_block;
}

std::shared_ptr<Mesh>
DistributedMeshGenerator::SetupLocalMesh(DistributedMeshData& mesh_info)
{
  auto grid_ptr = Mesh::New();
  for (auto& [id, name] : mesh_info.boundary_id_map)
    grid_ptr->SetBoundaryName(id, name);

  auto& vertices = mesh_info.vertices;
  for (const auto& [vid, vertex] : vertices)
    grid_ptr->AddGlobalVertex(vid, vertex);

  std::vector<Cell> local_cells;
  std::vector<Cell> ghost_cells;
  auto& cells = mesh_info.cells;
  for (auto& [pidgid, cell] : cells)
  {
    const auto& [cell_pid, cell_global_id] = pidgid;
    cell.global_id = cell_global_id;
    cell.partition_id = cell_pid;
    if (cell_pid == opensn::mpi_comm.rank())
      local_cells.push_back(std::move(cell));
    else
      ghost_cells.push_back(std::move(cell));
  }
  grid_ptr->SetCells(std::move(local_cells), std::move(ghost_cells), mesh_info.cell_connect);
  if (!mesh_info.cell_face_connect.empty())
    grid_ptr->SetCellFaces(mesh_info.cell_faces, mesh_info.cell_face_connect);

  grid_ptr->SetDimension(mesh_info.dimension);
  grid_ptr->SetCoordinateSystem(mesh_info.coord_sys);
  grid_ptr->SetType(mesh_info.mesh_type);
  grid_ptr->SetExtruded(mesh_info.extruded);
  grid_ptr->SetOrthoAttributes(mesh_info.ortho_attributes);
  grid_ptr->SetGlobalVertexCount(mesh_info.num_global_vertices);
  grid_ptr->ComputeGeometricInfo();

  ComputeAndPrintStats(grid_ptr);

  return grid_ptr;
}

} // namespace opensn
