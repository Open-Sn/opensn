// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_generator/split_file_mesh_generator.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/data_types/byte_array.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include <filesystem>

namespace opensn
{

SplitFileMeshGenerator::SplitFileMeshGenerator(const InputParameters& params)
  : MeshGenerator(params),
    num_partitions_(params.GetParamValue<int>("num_partitions")),
    split_mesh_dir_path_(params.GetParamValue<std::string>("split_mesh_dir_path")),
    file_prefix_(params.GetParamValue<std::string>("file_prefix")),
    read_only_(params.GetParamValue<bool>("read_only")),
    verbosity_level_(params.GetParamValue<int>("verbosity_level"))
{
}

std::shared_ptr<MeshContinuum>
SplitFileMeshGenerator::Execute()
{
  const auto num_mpi = mpi_comm.size();
  const auto num_partitions = num_mpi == 1 ? num_partitions_ : num_mpi;

  if (mpi_comm.rank() == 0 and (not read_only_))
  {
    // Execute all input generators
    // Note these could be empty
    std::shared_ptr<UnpartitionedMesh> current_umesh = nullptr;
    for (const auto& mesh_generator_ptr : inputs_)
      current_umesh = mesh_generator_ptr->GenerateUnpartitionedMesh(current_umesh);

    // Generate final umesh
    current_umesh = GenerateUnpartitionedMesh(current_umesh);

    log.Log() << "Writing split-mesh with " << num_partitions << " parts";
    const auto cell_pids = PartitionMesh(*current_umesh, num_partitions);
    WriteSplitMesh(cell_pids, *current_umesh, num_partitions);
    log.Log() << "Split-mesh with " << num_partitions << " parts successfully created";
  } // if home location

  // Other locations wait here for files to be written
  mpi_comm.barrier();

  std::shared_ptr<MeshContinuum> grid_ptr;
  if (mpi_comm.size() == num_partitions)
  {
    log.Log() << "Reading split-mesh";
    auto mesh_info = ReadSplitMesh();

    grid_ptr = SetupLocalMesh(mesh_info);

    log.Log() << "Done reading split-mesh files";
  }
  else
  {
    log.Log0Warning() << "After creating a split-mesh with mpi-processes < "
                         "num_parts the program will now auto terminate. This is not an error "
                         "and is the default behavior for the SplitFileMeshGenerator.\n";
    opensn::mpi_comm.abort(EXIT_SUCCESS);
  }

  mpi_comm.barrier();
  return grid_ptr;
}

void
SplitFileMeshGenerator::WriteSplitMesh(const std::vector<int>& cell_pids,
                                       const UnpartitionedMesh& umesh,
                                       const int num_partitions) const
{
  const std::filesystem::path dir_path = std::filesystem::absolute(split_mesh_dir_path_);

  const auto parent_path = dir_path.parent_path();
  if (not std::filesystem::exists(parent_path))
    throw std::runtime_error("Parent path " + parent_path.string() + " does not exist");

  bool root_dir_created = true;
  if (not std::filesystem::exists(dir_path))
    root_dir_created = std::filesystem::create_directories(dir_path);
  if (not root_dir_created)
    throw std::runtime_error("Failed to create directory " + dir_path.string() + ".");

  const auto& vertex_subs = umesh.GetVertextCellSubscriptions();
  const auto& raw_cells = umesh.GetRawCells();
  const auto& raw_vertices = umesh.GetVertices();

  uint64_t aux_counter = 0;
  for (int pid = 0; pid < num_partitions; ++pid)
  {
    const std::filesystem::path file_path =
      dir_path.string() + "/" + file_prefix_ + "_" + std::to_string(pid) + ".cmesh";

    std::ofstream ofile(file_path.string(), std::ios_base::binary | std::ios_base::out);
    if (not ofile.is_open())
      throw std::runtime_error("Failed to open " + file_path.string() + ".");

    // Appropriate cells and vertices to the current part being writting
    std::vector<uint64_t> local_cells_needed;
    std::set<uint64_t> cells_needed;
    std::set<uint64_t> vertices_needed;
    {
      local_cells_needed.reserve(raw_cells.size() / num_partitions);
      {
        uint64_t cell_global_id = 0;
        for (const auto cell_pid : cell_pids)
        {
          if (cell_pid == pid)
            local_cells_needed.push_back(cell_global_id);
          ++cell_global_id;
        }
      }

      for (const auto cell_global_id : local_cells_needed)
      {
        cells_needed.insert(cell_global_id);
        const auto& raw_cell = *raw_cells[cell_global_id];
        for (const auto vid : raw_cell.vertex_ids)
        {
          vertices_needed.insert(vid);
          for (const auto ghost_gid : vertex_subs[vid])
          {
            if (ghost_gid == cell_global_id)
              continue;

            cells_needed.insert(ghost_gid);
            const auto& ghost_raw_cell = *raw_cells[ghost_gid];
            for (const auto gvid : ghost_raw_cell.vertex_ids)
              vertices_needed.insert(gvid);
          }
        }
      }
    }

    if (verbosity_level_ >= 2)
      log.Log() << "Writing part " << pid << " num_local_cells=" << local_cells_needed.size();

    // Write mesh attributes and general info
    WriteBinaryValue(ofile, num_partitions); // int
    WriteBinaryValue<unsigned int>(ofile, umesh.GetDimension());
    WriteBinaryValue(ofile, static_cast<int>(umesh.GetCoordinateSystem()));
    WriteBinaryValue(ofile, static_cast<int>(umesh.GetType()));
    WriteBinaryValue(ofile, umesh.IsExtruded());
    const auto& [Nx, Ny, Nz] = umesh.GetOrthoAttributes();
    WriteBinaryValue(ofile, Nx); // size_t
    WriteBinaryValue(ofile, Ny); // size_t
    WriteBinaryValue(ofile, Nz); // size_t

    WriteBinaryValue(ofile, raw_vertices.size()); // size_t

    // Write the boundary map
    const auto& bndry_map = umesh.GetBoundaryIDMap();
    WriteBinaryValue(ofile, bndry_map.size()); // size_t
    for (const auto& [bid, bname] : bndry_map)
    {
      WriteBinaryValue(ofile, bid); // uint64_t
      const size_t num_chars = bname.size();
      WriteBinaryValue(ofile, num_chars);                     // size_t
      ofile.write(bname.data(), static_cast<int>(num_chars)); // characters
    }

    // Write how many cells and vertices in file
    WriteBinaryValue(ofile, cells_needed.size());    // size_t
    WriteBinaryValue(ofile, vertices_needed.size()); // size_t

    // Write cells
    constexpr size_t BUFFER_SIZE = 4096 * 2;
    ByteArray serial_data;
    serial_data.Data().reserve(BUFFER_SIZE * 2);
    for (const auto& cell_global_id : cells_needed)
    {
      const auto& cell = *raw_cells[cell_global_id];
      serial_data.Write(cell_pids[cell_global_id]);
      serial_data.Write(cell_global_id);
      SerializeCell(cell, serial_data);
      if (serial_data.Size() > BUFFER_SIZE)
      {
        WriteBinaryValue(ofile, serial_data);
        serial_data.Clear();
        serial_data.Data().reserve(serial_data.Data().capacity());
      }
    }
    if (serial_data.Size() > 0)
    {
      WriteBinaryValue(ofile, serial_data);

      serial_data.Clear();
    }

    // Write vertices
    for (const uint64_t vid : vertices_needed)
    {
      serial_data.Write(vid); // uint64_t
      serial_data.Write(raw_vertices[vid]);
      if (serial_data.Size() > BUFFER_SIZE)
      {
        WriteBinaryValue(ofile, serial_data);

        serial_data.Clear();
      }
    }
    if (serial_data.Size() > 0)
    {
      WriteBinaryValue(ofile, serial_data);

      serial_data.Clear();
    }
    ofile.close();

    const auto fraction_complete = static_cast<double>(pid) / static_cast<double>(num_partitions);
    if (fraction_complete >= static_cast<double>(aux_counter + 1) * 0.1)
    {
      if (verbosity_level_ >= 1)
        log.Log() << program_timer.GetTimeString() << " Surpassing part " << pid << " of "
                  << num_partitions << " (" << (aux_counter + 1) * 10 << "%)";
      ++aux_counter;
    }
  } // for p
}

SplitFileMeshGenerator::SplitMeshInfo
SplitFileMeshGenerator::ReadSplitMesh() const
{
  const auto pid = mpi_comm.rank();
  const std::filesystem::path dir_path = std::filesystem::absolute(split_mesh_dir_path_);
  const std::filesystem::path file_path =
    dir_path.string() + "/" + file_prefix_ + "_" + std::to_string(pid) + ".cmesh";

  SplitMeshInfo info_block;
  auto& cells = info_block.cells;
  auto& vertices = info_block.vertices;
  std::ifstream ifile(file_path, std::ios_base::binary | std::ios_base::in);
  if (not ifile.is_open())
    throw std::runtime_error("Failed to open " + file_path.string() + ".");

  // Read mesh attributes and general info
  const auto file_num_parts = ReadBinaryValue<int>(ifile);
  if (mpi_comm.size() != file_num_parts)
    throw std::logic_error("Split mesh files with prefix \"" + file_prefix_ +
                           "\" has been created with " + std::to_string(file_num_parts) +
                           " parts but is now being read with " + std::to_string(mpi_comm.size()) +
                           " processes.");

  info_block.dimension = ReadBinaryValue<unsigned int>(ifile);
  info_block.coord_sys = static_cast<CoordinateSystemType>(ReadBinaryValue<int>(ifile));
  info_block.mesh_type = static_cast<MeshType>(ReadBinaryValue<int>(ifile));
  info_block.extruded = ReadBinaryValue<bool>(ifile);
  info_block.ortho_attributes.Nx = ReadBinaryValue<size_t>(ifile);
  info_block.ortho_attributes.Ny = ReadBinaryValue<size_t>(ifile);
  info_block.ortho_attributes.Nz = ReadBinaryValue<size_t>(ifile);
  info_block.num_global_vertices = ReadBinaryValue<size_t>(ifile);

  // Read boundary map
  const auto num_boundaries = ReadBinaryValue<size_t>(ifile);
  for (size_t b = 0; b < num_boundaries; ++b)
  {
    const auto bid = ReadBinaryValue<uint64_t>(ifile);
    const auto num_chars = ReadBinaryValue<size_t>(ifile);
    std::string bname(num_chars, ' ');
    ifile.read(bname.data(), static_cast<int>(num_chars));
    info_block.boundary_id_map.insert(std::make_pair(bid, bname));
  }

  // Write how many cells and vertices in file
  const auto num_cells = ReadBinaryValue<size_t>(ifile);
  const auto num_vertices = ReadBinaryValue<size_t>(ifile);

  // Read the cells
  for (size_t c = 0; c < num_cells; ++c)
  {
    const auto cell_pid = ReadBinaryValue<int>(ifile);
    const auto cell_gid = ReadBinaryValue<uint64_t>(ifile);
    const auto cell_type = ReadBinaryValue<CellType>(ifile);
    const auto cell_sub_type = ReadBinaryValue<CellType>(ifile);

    UnpartitionedMesh::LightWeightCell new_cell(cell_type, cell_sub_type);

    new_cell.centroid = ReadBinaryValue<Vector3>(ifile);
    new_cell.block_id = ReadBinaryValue<int>(ifile);

    const auto num_vids = ReadBinaryValue<size_t>(ifile);
    for (size_t v = 0; v < num_vids; ++v)
      new_cell.vertex_ids.push_back(ReadBinaryValue<uint64_t>(ifile));

    const auto num_faces = ReadBinaryValue<size_t>(ifile);
    for (size_t f = 0; f < num_faces; ++f)
    {
      UnpartitionedMesh::LightWeightFace new_face;
      const auto num_face_vids = ReadBinaryValue<size_t>(ifile);
      for (size_t v = 0; v < num_face_vids; ++v)
        new_face.vertex_ids.push_back(ReadBinaryValue<uint64_t>(ifile));

      new_face.has_neighbor = ReadBinaryValue<bool>(ifile);
      new_face.neighbor = ReadBinaryValue<uint64_t>(ifile);

      new_cell.faces.push_back(std::move(new_face));
    } // for f

    cells.insert(std::make_pair(std::make_pair(cell_pid, cell_gid), std::move(new_cell)));
  } // for cell c

  // Read the vertices
  for (size_t v = 0; v < num_vertices; ++v)
  {
    const auto vid = ReadBinaryValue<uint64_t>(ifile);
    const auto vertex = ReadBinaryValue<Vector3>(ifile);
    vertices.insert(std::make_pair(vid, vertex));
  } // for vertex v

  ifile.close();
  return info_block;
}

OpenSnRegisterObjectInNamespace(mesh, SplitFileMeshGenerator);

InputParameters
SplitFileMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription(
    "Generates the mesh only on location 0, thereafter partitions the mesh"
    " but instead of broadcasting the mesh to other locations it creates binary"
    " mesh files for each location.");

  params.AddOptionalParameter("num_partitions",
                              0,
                              "The number of partitions to generate. If zero will "
                              "default to the number of MPI processes. Is "
                              "ignored if the number of MPI processes > 1.");

  params.AddOptionalParameter(
    "split_mesh_dir_path",
    "split_mesh",
    "Path of the directory to be created for containing the split meshes.");

  params.AddOptionalParameter(
    "file_prefix", input_path.stem().string(), "Prefix to use for all split mesh files");

  params.AddOptionalParameter(
    "read_only", false, "Controls whether the split mesh is recreated or just read.");

  params.AddOptionalParameter(
    "verbosity_level",
    1,
    "Verbosity level. 1 will report each 10% complete. 2 will print each part "
    "and the number of local cells it wrote.");

  return params;
}

std::shared_ptr<SplitFileMeshGenerator>
SplitFileMeshGenerator::Create(const ParameterBlock& params)
{
  auto& factory = ObjectFactory::GetInstance();
  return factory.Create<SplitFileMeshGenerator>("mesh::SplitFileMeshGenerator", params);
}

std::shared_ptr<MeshContinuum>
SplitFileMeshGenerator::SetupLocalMesh(SplitMeshInfo& mesh_info)
{
  auto grid_ptr = MeshContinuum::New();
  grid_ptr->GetBoundaryIDMap() = mesh_info.boundary_id_map;

  auto& cells = mesh_info.cells;
  auto& vertices = mesh_info.vertices;

  for (const auto& [vid, vertex] : vertices)
    grid_ptr->vertices.Insert(vid, vertex);

  for (const auto& [pidgid, raw_cell] : cells)
  {
    const auto& [cell_pid, cell_global_id] = pidgid;
    auto cell = SetupCell(raw_cell, cell_global_id, cell_pid);
    grid_ptr->cells.PushBack(std::move(cell));
  }

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

void
SplitFileMeshGenerator::SerializeCell(const UnpartitionedMesh::LightWeightCell& cell,
                                      ByteArray& serial_buffer)
{
  serial_buffer.Write(cell.type);
  serial_buffer.Write(cell.sub_type);
  serial_buffer.Write(cell.centroid);
  serial_buffer.Write(cell.block_id);
  serial_buffer.Write(cell.vertex_ids.size());
  for (const uint64_t vid : cell.vertex_ids)
    serial_buffer.Write(vid);
  serial_buffer.Write(cell.faces.size());
  for (const auto& face : cell.faces)
  {
    serial_buffer.Write(face.vertex_ids.size());
    for (const uint64_t vid : face.vertex_ids)
      serial_buffer.Write(vid);
    serial_buffer.Write(face.has_neighbor);
    serial_buffer.Write(face.neighbor);
  }
}

} // namespace opensn
