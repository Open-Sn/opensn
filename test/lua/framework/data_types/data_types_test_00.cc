#include "lua/lib/console.h"
#include "framework/data_types/byte_array.h"
#include "framework/data_types/ndarray.h"
#include "framework/mesh/cell/cell.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi_utils.h"
#include <map>

using namespace opensn;

namespace unit_tests
{

void
data_types_Test00()
{
  bool passed = true;

  opensn::log.Log() << "GOLD_BEGIN";

  // Testing Byte array
  // serialization
  opensn::log.Log() << "Testing ByteArray "
                       "Serialization/DeSerialization\n";
  if (opensn::mpi_comm.size() == 2)
  {

    std::map<int /*pid*/, ByteArray> send_data;

    Cell poster_child_cell(CellType::POLYHEDRON, CellType::HEXAHEDRON);
    {
      poster_child_cell.global_id = 321;
      poster_child_cell.local_id = 123;
      poster_child_cell.partition_id = 0;
      poster_child_cell.centroid = Vector3(0.5, 0.5, 0.5);
      poster_child_cell.block_id = -2;

      poster_child_cell.vertex_ids = {0, 1, 2, 3, 4, 5, 6, 7};

      // Bottom face
      {
        CellFace face;
        face.vertex_ids = {0, 3, 2, 1};
        face.normal = {0, 0, -1};
        face.centroid = {0.5, 0.5, 0.0};
        face.has_neighbor = false;
        face.neighbor_id = 0;
        poster_child_cell.faces.push_back(std::move(face));
      }
      // Top face
      {
        CellFace face;
        face.vertex_ids = {4, 5, 6, 7};
        face.normal = {0, 0, 1};
        face.centroid = {0.5, 0.5, 1.0};
        face.has_neighbor = false;
        face.neighbor_id = 1;
        poster_child_cell.faces.push_back(std::move(face));
      }
      // Left face
      {
        CellFace face;
        face.vertex_ids = {0, 4, 7, 3};
        face.normal = {-1, 0, 0};
        face.centroid = {0.0, 0.5, 0.5};
        face.has_neighbor = false;
        face.neighbor_id = 2;
        poster_child_cell.faces.push_back(std::move(face));
      }
      // Right face
      {
        CellFace face;
        face.vertex_ids = {1, 2, 6, 5};
        face.normal = {1, 0, 0};
        face.centroid = {1.0, 0.5, 0.5};
        face.has_neighbor = false;
        face.neighbor_id = 3;
        poster_child_cell.faces.push_back(std::move(face));
      }
      // Front face
      {
        CellFace face;
        face.vertex_ids = {0, 1, 5, 4};
        face.normal = {0, -1, 0};
        face.centroid = {0.5, 0.0, 0.5};
        face.has_neighbor = false;
        face.neighbor_id = 4;
        poster_child_cell.faces.push_back(std::move(face));
      }
      // Back face
      {
        CellFace face;
        face.vertex_ids = {3, 7, 6, 2};
        face.normal = {0, 1, 0};
        face.centroid = {0.5, 1.0, 0.5};
        face.has_neighbor = false;
        face.neighbor_id = 5;
        poster_child_cell.faces.push_back(std::move(face));
      }
    }

    if (opensn::mpi_comm.rank() == 0)
    {
      send_data[1].Append(poster_child_cell.Serialize());
      send_data[1].Append(poster_child_cell.Serialize().Data());
    }

    std::map<int /*pid*/, std::vector<std::byte>> send_data_bytes;

    for (const auto& pid_byte_array : send_data)
      send_data_bytes[pid_byte_array.first] = pid_byte_array.second.Data();

    std::map<int /*pid*/, std::vector<std::byte>> recv_data_bytes = MapAllToAll(send_data_bytes);

    for (const auto& pid_vec_bytes : recv_data_bytes)
    {
      // auto& pid = pid_vec_bytes.first;
      auto& vec_bytes = pid_vec_bytes.second;

      ByteArray byte_array(vec_bytes);

      size_t address = 0;
      while (address < byte_array.Size())
      {
        const Cell read_cell = Cell::DeSerialize(byte_array, address);

        auto& rcell = read_cell;
        auto& pcell = poster_child_cell;

        if (rcell.GetType() != pcell.GetType())
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.GetSubType() != pcell.GetSubType())
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.global_id != pcell.global_id)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.local_id != pcell.local_id)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.partition_id != pcell.partition_id)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.block_id != pcell.block_id)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.vertex_ids != pcell.vertex_ids)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }

        if (rcell.faces.size() != pcell.faces.size())
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }

        size_t f = 0;
        for (const auto& rface : rcell.faces)
        {
          const auto& pface = pcell.faces[f];

          if (rface.vertex_ids != pface.vertex_ids)
          {
            passed = false;
            opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
            break;
          }
          if (rface.has_neighbor != pface.has_neighbor)
          {
            passed = false;
            opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
            break;
          }
          if (rface.neighbor_id != pface.neighbor_id)
          {
            passed = false;
            opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
            break;
          }
          ++f;
        }
      }
    }
  } // if #procs=2
  else
    throw std::invalid_argument("unit_tests::Test_data_types requires"
                                " at least 2 MPI processes.");

  if (not passed)
  {
    opensn::log.Log() << "ByteArray "
                         "Serialization/DeSerialization ... Failed\n";
  }
  else
    opensn::log.Log() << "ByteArray"
                         "Serialization/DeSerialization ... Passed\n";

  opensn::log.Log() << "GOLD_END";
}

BIND_FUNCTION(unit_tests, data_types_Test00);

} //  namespace unit_tests
