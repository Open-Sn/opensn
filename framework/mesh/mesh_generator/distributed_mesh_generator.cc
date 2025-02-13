// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_generator/distributed_mesh_generator.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/data_types/byte_array.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(mesh, DistributedMeshGenerator);

InputParameters
DistributedMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription(
    "Generates and partitions the mesh on location 0. The partitioned mesh is "
    "broadcast to all other locations.");
  params.SetDocGroup("doc_MeshGenerators");

  return params;
}

std::shared_ptr<DistributedMeshGenerator>
DistributedMeshGenerator::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<DistributedMeshGenerator>("mesh::DistributedMeshGenerator", params);
}

DistributedMeshGenerator::DistributedMeshGenerator(const InputParameters& params)
  : MeshGenerator(params), num_parts_(opensn::mpi_comm.size())
{
}

std::shared_ptr<MeshContinuum>
DistributedMeshGenerator::Execute()
{
  const int rank = opensn::mpi_comm.rank();
  const int num_parts = opensn::mpi_comm.size();
  DistributedMeshData mesh_info;

  log.Log() << program_timer.GetTimeString() << " Distributing mesh with " << num_parts << " parts";

  if (rank == 0)
  {
    std::shared_ptr<UnpartitionedMesh> current_umesh = nullptr;
    for (auto mesh_generator_ptr : inputs_)
    {
      auto new_umesh = mesh_generator_ptr->GenerateUnpartitionedMesh(current_umesh);
      current_umesh = new_umesh;
    }
    current_umesh = GenerateUnpartitionedMesh(current_umesh);

    const auto cell_pids = PartitionMesh(*current_umesh, num_parts);
    auto serial_data = DistributeSerializedMeshData(cell_pids, *current_umesh, num_parts);
    mesh_info = DeserializeMeshData(serial_data);
  }
  else
  {
    std::vector<std::byte> data;
    opensn::mpi_comm.recv<std::byte>(0, rank, data);
    ByteArray serial_data(data);
    mesh_info = DeserializeMeshData(serial_data);
  }

  auto grid_ptr = SetupLocalMesh(mesh_info);

  opensn::mpi_comm.barrier();

  log.Log() << program_timer.GetTimeString() << " Mesh successfully distributed";

  return grid_ptr;
}

ByteArray
DistributedMeshGenerator::DistributeSerializedMeshData(const std::vector<int64_t>& cell_pids,
                                                       const UnpartitionedMesh& umesh,
                                                       int num_parts)
{
  const auto& vertex_subs = umesh.GetVertextCellSubscriptions();
  const auto& raw_cells = umesh.GetRawCells();
  const auto& raw_vertices = umesh.GetVertices();
  ByteArray loc0_data;

  for (int pid = 0; pid < num_parts; ++pid)
  {
    ByteArray serial_data;

    std::vector<uint64_t> local_cells_needed;
    std::set<uint64_t> cells_needed;
    std::set<uint64_t> vertices_needed;
    local_cells_needed.reserve(cell_pids.size() / num_parts);

    for (uint64_t cell_global_id = 0; cell_global_id < cell_pids.size(); ++cell_global_id)
    {
      if (cell_pids[cell_global_id] == pid)
      {
        const auto& raw_cell = *raw_cells[cell_global_id];
        local_cells_needed.push_back(cell_global_id);
        cells_needed.emplace(cell_global_id);

        for (uint64_t vid : raw_cell.vertex_ids)
        {
          vertices_needed.emplace(vid);

          // Process ghost cells
          for (uint64_t ghost_gid : vertex_subs[vid])
          {
            if (ghost_gid != cell_global_id && cells_needed.find(ghost_gid) == cells_needed.end())
            {
              cells_needed.emplace(ghost_gid);
              const auto& ghost_raw_cell = *raw_cells[ghost_gid];
              // Insert ghost vertex IDs
              for (uint64_t gvid : ghost_raw_cell.vertex_ids)
                vertices_needed.emplace(gvid);
            }
          }
        }
      }
    }

    // Basic mesh data
    serial_data.Write<unsigned int>(umesh.GetDimension());
    serial_data.Write(static_cast<int>(umesh.GetType()));
    serial_data.Write(umesh.IsExtruded());
    auto& ortho_attrs = umesh.GetOrthoAttributes();
    serial_data.Write(ortho_attrs.Nx);
    serial_data.Write(ortho_attrs.Ny);
    serial_data.Write(ortho_attrs.Nz);
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
        serial_data.Write(bname.data()[i]);
    }

    // Number of cells and vertices
    serial_data.Write(cells_needed.size());
    serial_data.Write(vertices_needed.size());

    // Cell data
    for (const auto& cell_global_id : cells_needed)
    {
      const auto& cell = *raw_cells[cell_global_id];
      serial_data.Write(static_cast<int>(cell_pids[cell_global_id]));
      serial_data.Write(cell_global_id);
      serial_data.Write(cell.type);
      serial_data.Write(cell.sub_type);
      serial_data.Write(cell.centroid.x);
      serial_data.Write(cell.centroid.y);
      serial_data.Write(cell.centroid.z);
      serial_data.Write(cell.material_id);
      serial_data.Write(cell.vertex_ids.size());
      for (uint64_t vid : cell.vertex_ids)
        serial_data.Write(vid);

      serial_data.Write(cell.faces.size());
      for (const auto& face : cell.faces)
      {
        serial_data.Write(face.vertex_ids.size());
        for (uint64_t vid : face.vertex_ids)
          serial_data.Write(vid);
        serial_data.Write(face.has_neighbor);
        serial_data.Write(face.neighbor);
      }
    }

    // Vertex data
    for (uint64_t vid : vertices_needed)
    {
      serial_data.Write(vid);
      serial_data.Write(raw_vertices[vid].x);
      serial_data.Write(raw_vertices[vid].y);
      serial_data.Write(raw_vertices[vid].z);
    }

    if (pid == 0)
      loc0_data = serial_data;
    else
      opensn::mpi_comm.send<std::byte>(pid, pid, serial_data.Data().data(), serial_data.Size());
  }

  return loc0_data;
}

DistributedMeshGenerator::DistributedMeshData
DistributedMeshGenerator::DeserializeMeshData(ByteArray& serial_data)
{
  DistributedMeshData info_block;

  // Basic mesh data
  info_block.dimension = serial_data.Read<unsigned int>();
  info_block.mesh_type = static_cast<MeshType>(serial_data.Read<int>());
  info_block.extruded = serial_data.Read<bool>();
  info_block.ortho_attributes.Nx = serial_data.Read<size_t>();
  info_block.ortho_attributes.Ny = serial_data.Read<size_t>();
  info_block.ortho_attributes.Nz = serial_data.Read<size_t>();
  info_block.num_global_vertices = serial_data.Read<size_t>();

  // Boundaries
  auto num_bndries = serial_data.Read<size_t>();
  for (auto b = 0; b < num_bndries; ++b)
  {
    auto bid = serial_data.Read<uint64_t>();
    auto num_chars = serial_data.Read<size_t>();
    std::string bname(num_chars, ' ');
    for (auto i = 0; i < num_chars; ++i)
      bname.data()[i] = serial_data.Read<char>();
    info_block.boundary_id_map.insert(std::make_pair(bid, bname));
  }

  // Number of cells and vertices
  auto num_cells = serial_data.Read<size_t>();
  auto num_vertices = serial_data.Read<size_t>();

  // Cell data
  for (size_t i = 0; i < num_cells; ++i)
  {
    int cell_pid = serial_data.Read<int>();
    auto cell_gid = serial_data.Read<uint64_t>();
    auto type = serial_data.Read<CellType>();
    auto sub_type = serial_data.Read<CellType>();

    UnpartitionedMesh::LightWeightCell cell(type, sub_type);

    cell.centroid.x = serial_data.Read<double>();
    cell.centroid.y = serial_data.Read<double>();
    cell.centroid.z = serial_data.Read<double>();
    cell.material_id = serial_data.Read<int>();

    auto num_vids = serial_data.Read<size_t>();
    for (size_t v = 0; v < num_vids; ++v)
      cell.vertex_ids.push_back(serial_data.Read<uint64_t>());

    auto num_faces = serial_data.Read<size_t>();
    for (size_t f = 0; f < num_faces; ++f)
    {
      UnpartitionedMesh::LightWeightFace face;
      auto num_face_vids = serial_data.Read<size_t>();
      for (size_t v = 0; v < num_face_vids; ++v)
        face.vertex_ids.push_back(serial_data.Read<uint64_t>());

      face.has_neighbor = serial_data.Read<bool>();
      face.neighbor = serial_data.Read<uint64_t>();

      cell.faces.push_back(std::move(face));
    }
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

std::shared_ptr<MeshContinuum>
DistributedMeshGenerator::SetupLocalMesh(DistributedMeshData& mesh_info)
{
  auto grid_ptr = MeshContinuum::New();
  grid_ptr->GetBoundaryIDMap() = mesh_info.boundary_id_map;

  auto& vertices = mesh_info.vertices;
  for (const auto& [vid, vertex] : vertices)
    grid_ptr->vertices.Insert(vid, vertex);

  auto& cells = mesh_info.cells;
  for (const auto& [pidgid, raw_cell] : cells)
  {
    const auto& [cell_pid, cell_global_id] = pidgid;
    auto cell = SetupCell(raw_cell, cell_global_id, cell_pid, STLVertexListHelper(vertices));

    grid_ptr->cells.PushBack(std::move(cell));
  }

  grid_ptr->SetDimension(mesh_info.dimension);
  grid_ptr->SetType(mesh_info.mesh_type);
  grid_ptr->SetExtruded(mesh_info.extruded);
  grid_ptr->SetOrthoAttributes(mesh_info.ortho_attributes);
  grid_ptr->SetGlobalVertexCount(mesh_info.num_global_vertices);
  ComputeAndPrintStats(*grid_ptr);

  return grid_ptr;
}

} // namespace opensn
