#include "framework/data_types/byte_array.h"
#include "framework/data_types/ndarray.h"
#include "framework/mesh/Cell/cell.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "framework/console/chi_console.h"

#include "framework/mpi/chi_mpi.h"

#include "framework/mpi/chi_mpi_utils_map_all2all.h"

#include <map>

namespace chi_unit_tests
{

chi::ParameterBlock chi_data_types_Test00(const chi::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chi_data_types_Test00,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/chi_data_types_Test00);

chi::ParameterBlock
chi_data_types_Test00(const chi::InputParameters&)
{
  bool passed = true;

  //======================================================= Byte array
  // write/read
  Chi::log.Log() << "GOLD_BEGIN";
  Chi::log.Log() << "Testing chi_data_types::ByteArray Write and Read\n";
  chi_data_types::ByteArray barr;

  barr.Write<double>(1.01234567890123456789);
  barr.Write<int>(-15600700);
  barr.Write<double>(2.0987654321);
  barr.Write<bool>(false);
  barr.Write<bool>(true);
  barr.Write<chi_mesh::Vector3>(chi_mesh::Vector3(3.0, 2.0, 1.0));

  auto dbl_value = barr.Read<double>();
  auto int_value = barr.Read<int>();
  auto dbl_value2 = barr.Read<double>();
  auto bl_value1 = barr.Read<bool>();
  auto bl_value2 = barr.Read<bool>();
  auto vec3 = barr.Read<chi_mesh::Vector3>();

  chi_data_types::ByteArray seeker;
  seeker.Write<bool>(false);
  seeker.Write<double>(1.01234567890123456789);

  Chi::log.Log() << "EndOfBuffer " << seeker.EndOfBuffer();
  Chi::log.Log() << "Offset " << seeker.Offset();
  seeker.Seek(seeker.Size() - sizeof(double));
  Chi::log.Log() << "OffsetAfterSeek " << seeker.Offset();
  Chi::log.Log() << "Value check " << seeker.Read<double>();

  if (dbl_value != 1.01234567890123456789 or int_value != -15600700 or dbl_value2 != 2.0987654321 or
      bl_value1 or !bl_value2 or vec3[0] != chi_mesh::Vector3(3.0, 2.0, 1.0)[0] or
      vec3[1] != chi_mesh::Vector3(3.0, 2.0, 1.0)[1] or
      vec3[2] != chi_mesh::Vector3(3.0, 2.0, 1.0)[2])

  {
    passed = false;
    Chi::log.Log() << std::string("chi_data_types::ByteArray"
                                  " Write/Read ... Failed\n");
  }
  else
    Chi::log.Log() << std::string("chi_data_types::ByteArray "
                                  "Write/Read ... Passed\n");

  //======================================================= Testing Byte array
  // serialization
  Chi::log.Log() << "Testing chi_data_types::ByteArray "
                    "Serialization/DeSerialization\n";
  if (Chi::mpi.process_count == 2)
  {

    std::map<int /*pid*/, chi_data_types::ByteArray> send_data;

    chi_mesh::Cell poster_child_cell(chi_mesh::CellType::POLYHEDRON,
                                     chi_mesh::CellType::HEXAHEDRON);
    {
      poster_child_cell.global_id_ = 321;
      poster_child_cell.local_id_ = 123;
      poster_child_cell.partition_id_ = 0;
      poster_child_cell.centroid_ = chi_mesh::Vector3(0.5, 0.5, 0.5);
      poster_child_cell.material_id_ = -2;

      poster_child_cell.vertex_ids_ = {0, 1, 2, 3, 4, 5, 6, 7};

      // Bottom face
      {
        chi_mesh::CellFace face;
        face.vertex_ids_ = {0, 3, 2, 1};
        face.normal_ = {0, 0, -1};
        face.centroid_ = {0.5, 0.5, 0.0};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 0;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Top face
      {
        chi_mesh::CellFace face;
        face.vertex_ids_ = {4, 5, 6, 7};
        face.normal_ = {0, 0, 1};
        face.centroid_ = {0.5, 0.5, 1.0};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 1;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Left face
      {
        chi_mesh::CellFace face;
        face.vertex_ids_ = {0, 4, 7, 3};
        face.normal_ = {-1, 0, 0};
        face.centroid_ = {0.0, 0.5, 0.5};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 2;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Right face
      {
        chi_mesh::CellFace face;
        face.vertex_ids_ = {1, 2, 6, 5};
        face.normal_ = {1, 0, 0};
        face.centroid_ = {1.0, 0.5, 0.5};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 3;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Front face
      {
        chi_mesh::CellFace face;
        face.vertex_ids_ = {0, 1, 5, 4};
        face.normal_ = {0, -1, 0};
        face.centroid_ = {0.5, 0.0, 0.5};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 4;
        poster_child_cell.faces_.push_back(std::move(face));
      }
      // Back face
      {
        chi_mesh::CellFace face;
        face.vertex_ids_ = {3, 7, 6, 2};
        face.normal_ = {0, 1, 0};
        face.centroid_ = {0.5, 1.0, 0.5};
        face.has_neighbor_ = false;
        face.neighbor_id_ = 5;
        poster_child_cell.faces_.push_back(std::move(face));
      }
    }

    if (Chi::mpi.location_id == 0)
    {
      send_data[1].Append(poster_child_cell.Serialize());
      send_data[1].Append(poster_child_cell.Serialize().Data());
    }

    std::map<int /*pid*/, std::vector<std::byte>> send_data_bytes;

    for (const auto& pid_byte_array : send_data)
      send_data_bytes[pid_byte_array.first] = pid_byte_array.second.Data();

    std::map<int /*pid*/, std::vector<std::byte>> recv_data_bytes =
      chi_mpi_utils::MapAllToAll(send_data_bytes, MPI_BYTE);

    for (const auto& pid_vec_bytes : recv_data_bytes)
    {
      // auto& pid = pid_vec_bytes.first;
      auto& vec_bytes = pid_vec_bytes.second;

      chi_data_types::ByteArray byte_array(vec_bytes);

      size_t address = 0;
      while (address < byte_array.Size())
      {
        const chi_mesh::Cell read_cell = chi_mesh::Cell::DeSerialize(byte_array, address);

        auto& rcell = read_cell;
        auto& pcell = poster_child_cell;

        if (rcell.Type() != pcell.Type())
        {
          passed = false;
          Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.SubType() != pcell.SubType())
        {
          passed = false;
          Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.global_id_ != pcell.global_id_)
        {
          passed = false;
          Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.local_id_ != pcell.local_id_)
        {
          passed = false;
          Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.partition_id_ != pcell.partition_id_)
        {
          passed = false;
          Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.material_id_ != pcell.material_id_)
        {
          passed = false;
          Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }
        if (rcell.vertex_ids_ != pcell.vertex_ids_)
        {
          passed = false;
          Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }

        if (rcell.faces_.size() != pcell.faces_.size())
        {
          passed = false;
          Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
          break;
        }

        size_t f = 0;
        for (const auto& rface : rcell.faces_)
        {
          const auto& pface = pcell.faces_[f];

          if (rface.vertex_ids_ != pface.vertex_ids_)
          {
            passed = false;
            Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
            break;
          }
          if (rface.has_neighbor_ != pface.has_neighbor_)
          {
            passed = false;
            Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
            break;
          }
          if (rface.neighbor_id_ != pface.neighbor_id_)
          {
            passed = false;
            Chi::log.Log0Error() << "Line: " << __LINE__ << "\n";
            break;
          }
          ++f;
        }
      }
    }
  } // if #procs=2
  else
    throw std::invalid_argument("chi_unit_tests::Test_chi_data_types requires"
                                " at least 2 MPI processes.");

  if (not passed)
  {
    Chi::log.Log() << "chi_data_types::ByteArray "
                      "Serialization/DeSerialization ... Failed\n";
  }
  else
    Chi::log.Log() << "chi_data_types::ByteArray"
                      "Serialization/DeSerialization ... Passed\n";

  //======================================================= Testing NDArray
  //
  Chi::log.Log() << "Testing chi_data_types::NDArray\n";
  std::stringstream dummy;
  // Constructor vector
  // rank()
  {
    chi_data_types::NDArray<double> nd_array1(std::vector<size_t>{2, 2, 2});
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
    chi_data_types::NDArray<double> nd_array1(std::array<size_t, 3>{2, 2, 2});
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
    chi_data_types::NDArray<double> nd_array1({2, 2, 2});
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done3\n";
  }

  // Constructor vector value
  {
    dummy << "Should be 2x2x2=8 zeros\n";
    chi_data_types::NDArray<double> nd_array1(std::vector<size_t>{2, 2, 2}, 0.0);
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done4\n";
  }

  // Constructor array value
  {
    dummy << "Should be 2x2x2=8 zeros\n";
    chi_data_types::NDArray<double> nd_array1(std::array<size_t, 3>{2, 2, 2}, 0.0);
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done5\n";
  }

  // Constructor initializer_list value
  {
    dummy << "Should be 2x2x2=8 zeros\n";
    chi_data_types::NDArray<double> nd_array1({2, 2, 2}, 0.0);
    for (auto val : nd_array1)
    {
      dummy << val << " ";
    }
    dummy << "Done6\n";
  }

  // Constructor none
  {
    dummy << "Should not print anything\n";
    chi_data_types::NDArray<double> nd_array1;
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
    chi_data_types::NDArray<double> nd_array2(std::vector<size_t>{2, 2, 2});
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
    chi_data_types::NDArray<double> nd_array4(std::array<size_t, 3>{2, 2, 2});
    nd_array4.set(1.0);

    dummy << "size " << nd_array4.size() << "\n";
    dummy << "empty() " << nd_array4.empty();

    dummy << "Should be 2x2x2=8 ones\n";
    for (auto val : nd_array4)
      dummy << val << " ";
    dummy << "\n";
    dummy << "Done10\n";
  }

  Chi::log.Log() << dummy.str();

  Chi::log.Log() << "GOLD_END";

  return chi::ParameterBlock();
}

} // namespace chi_unit_tests
