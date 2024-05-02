#include "framework/data_types/byte_array.h"
#include "framework/data_types/ndarray.h"
#include "framework/mesh/cell/cell.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "lua/framework/console/console.h"

#include "framework/mpi/mpi_utils.h"

#include <map>

using namespace opensn;

namespace unit_tests
{

ParameterBlock data_types_Test00(const InputParameters& params);

RegisterWrapperFunctionInNamespace(unit_tests, data_types_Test00, nullptr, data_types_Test00);

ParameterBlock
data_types_Test00(const InputParameters&)
{
  bool passed = true;

  // Byte array
  // write/read
  opensn::log.Log() << "GOLD_BEGIN";
  opensn::log.Log() << "Testing ByteArray Write and Read\n";
  ByteArray barr;

  barr.Write<double>(1.01234567890123456789);
  barr.Write<int>(-15600700);
  barr.Write<double>(2.0987654321);
  barr.Write<bool>(false);
  barr.Write<bool>(true);
  barr.Write<Vector3>(Vector3(3.0, 2.0, 1.0));

  auto dbl_value = barr.Read<double>();
  auto int_value = barr.Read<int>();
  auto dbl_value2 = barr.Read<double>();
  auto bl_value1 = barr.Read<bool>();
  auto bl_value2 = barr.Read<bool>();
  auto vec3 = barr.Read<Vector3>();

  ByteArray seeker;
  seeker.Write<bool>(false);
  seeker.Write<double>(1.01234567890123456789);

  opensn::log.Log() << "EndOfBuffer " << seeker.EndOfBuffer();
  opensn::log.Log() << "Offset " << seeker.Offset();
  seeker.Seek(seeker.Size() - sizeof(double));
  opensn::log.Log() << "OffsetAfterSeek " << seeker.Offset();
  opensn::log.Log() << "Value check " << seeker.Read<double>();

  if (dbl_value != 1.01234567890123456789 or int_value != -15600700 or dbl_value2 != 2.0987654321 or
      bl_value1 or !bl_value2 or vec3[0] != Vector3(3.0, 2.0, 1.0)[0] or
      vec3[1] != Vector3(3.0, 2.0, 1.0)[1] or vec3[2] != Vector3(3.0, 2.0, 1.0)[2])

  {
    passed = false;
    opensn::log.Log() << std::string("ByteArray"
                                     " Write/Read ... Failed\n");
  }
  else
    opensn::log.Log() << std::string("ByteArray "
                                     "Write/Read ... Passed\n");

  // Testing Byte array
  // serialization
  opensn::log.Log() << "Testing ByteArray "
                       "Serialization/DeSerialization\n";
  if (opensn::mpi_comm.size() == 2)
  {

    std::map<int /*pid*/, ByteArray> send_data;

    Cell poster_child_cell(CellType::POLYHEDRON, CellType::HEXAHEDRON);
    {
      poster_child_cell.global_id_ = 321;
      poster_child_cell.local_id_ = 123;
      poster_child_cell.partition_id_ = 0;
      poster_child_cell.centroid_ = Vector3(0.5, 0.5, 0.5);
      poster_child_cell.material_id_ = -2;

      poster_child_cell.vertex_ids_ = {0, 1, 2, 3, 4, 5, 6, 7};

      // Bottom face
      {
        CellFace face;
        face.vertex_ids_ = {0, 3, 2, 1};
        face.normal_ = {0, 0, -1};
        face.centroid_ = {0.5, 0.5, 0.0};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 0;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Top face
      {
        CellFace face;
        face.vertex_ids_ = {4, 5, 6, 7};
        face.normal_ = {0, 0, 1};
        face.centroid_ = {0.5, 0.5, 1.0};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 1;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Left face
      {
        CellFace face;
        face.vertex_ids_ = {0, 4, 7, 3};
        face.normal_ = {-1, 0, 0};
        face.centroid_ = {0.0, 0.5, 0.5};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 2;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Right face
      {
        CellFace face;
        face.vertex_ids_ = {1, 2, 6, 5};
        face.normal_ = {1, 0, 0};
        face.centroid_ = {1.0, 0.5, 0.5};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 3;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Front face
      {
        CellFace face;
        face.vertex_ids_ = {0, 1, 5, 4};
        face.normal_ = {0, -1, 0};
        face.centroid_ = {0.5, 0.0, 0.5};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 4;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Back face
      {
        CellFace face;
        face.vertex_ids_ = {3, 7, 6, 2};
        face.normal_ = {0, 1, 0};
        face.centroid_ = {0.5, 1.0, 0.5};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 5;
        poster_child_cell.faces_.push_back(std::move(face));
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

        if (rcell.Type() != pcell.Type())
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.SubType() != pcell.SubType())
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.global_id_ != pcell.global_id_)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.local_id_ != pcell.local_id_)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.partition_id_ != pcell.partition_id_)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.material_id_ != pcell.material_id_)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.vertex_ids_ != pcell.vertex_ids_)
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }

        if (rcell.faces_.size() != pcell.faces_.size())
        {
          passed = false;
          opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }

        size_t f = 0;
        for (const auto& rface : rcell.faces_)
        {
          const auto& pface = pcell.faces_[f];

          if (rface.vertex_ids_ != pface.vertex_ids_)
          {
            passed = false;
            opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
            break;
          }
          if (rface.has_neighbor_ != pface.has_neighbor_)
          {
            passed = false;
            opensn::log.Log0Error() << "Line: " << __LINE__ << "\n";
            break;
          }
          if (rface.neighbor_id_ != pface.neighbor_id_)
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

  // Testing NDArray
  //
  opensn::log.Log() << "Testing NDArray\n";
  std::stringstream dummy;
  // Constructor vector
  // rank()
  {
    NDArray<double> nd_array1(std::vector<size_t>{2, 2, 2});
    dummy << "Should be printing rank and 2x2x2=8 zeros\n";
    dummy << nd_array1.rank() << "\n";
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done1\n";
  }

  // Constructor array
  {
    NDArray<double> nd_array1(std::array<size_t, 3>{2, 2, 2});
    dummy << "Should be 2x2x2=8 zeros\n";
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done2\n";
  }

  // Constructor initializer_list
  {
    dummy << "Should be 2x2x2=8 zeros\n";
    NDArray<double> nd_array1({2, 2, 2});
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done3\n";
  }

  // Constructor vector value
  {
    dummy << "Should be 2x2x2=8 zeros\n";
    NDArray<double> nd_array1(std::vector<size_t>{2, 2, 2}, 0.0);
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done4\n";
  }

  // Constructor array value
  {
    dummy << "Should be 2x2x2=8 zeros\n";
    NDArray<double> nd_array1(std::array<size_t, 3>{2, 2, 2}, 0.0);
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done5\n";
  }

  // Constructor initializer_list value
  {
    dummy << "Should be 2x2x2=8 zeros\n";
    NDArray<double> nd_array1({2, 2, 2}, 0.0);
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done6\n";
  }

  // Constructor none
  {
    dummy << "Should not print anything\n";
    NDArray<double> nd_array1;
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done7\n";
  }

  // method set
  // iterators
  // const iterators
  {
    NDArray<double> nd_array2(std::vector<size_t>{2, 2, 2});
    nd_array2.set(1.0);

    dummy << "Should be 2x2x2=8 ones\n";
    for (auto val : nd_array2)
      dummy << val << " ";
    dummy << "\n";
    dummy << "Done8\n";

    dummy << "Should be 2x2x2=8 ones\n";
    const auto& nd_array3 = nd_array2;
    for (auto i = nd_array3.cbegin(); i != nd_array3.cend(); ++i)
      dummy << *i << " ";
    dummy << "\n";
    dummy << "Done9\n";
  }

  // size and empty
  {
    NDArray<double> nd_array4(std::array<size_t, 3>{2, 2, 2});
    nd_array4.set(1.0);

    dummy << "size " << nd_array4.size() << "\n";
    dummy << "empty() " << nd_array4.empty();

    dummy << "Should be 2x2x2=8 ones\n";
    for (auto val : nd_array4)
      dummy << val << " ";
    dummy << "\n";
    dummy << "Done10\n";
  }

  opensn::log.Log() << dummy.str();

  opensn::log.Log() << "GOLD_END";

  return ParameterBlock();
}

} //  namespace unit_tests
