// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstdint>

#if defined(__NVCC__)
#define __inline_host_dev__ inline __host__ __device__
#else
#define __inline_host_dev__ inline
#endif

namespace opensn
{

/// Face view from contiguous block of memory
struct FaceView
{
  __inline_host_dev__ FaceView() {}
  __inline_host_dev__ void Update(const char* face_data)
  {
    // number of face nodes
    const std::uint64_t* num_face_nodes_data = reinterpret_cast<const std::uint64_t*>(face_data);
    num_face_nodes = *(num_face_nodes_data++);
    face_data = reinterpret_cast<const char*>(num_face_nodes_data);
    // outflow
    double* const* outflow_data = reinterpret_cast<double* const*>(face_data);
    outflow = *(outflow_data++);
    face_data = reinterpret_cast<const char*>(outflow_data);
    // normal
    const double* normal_vector_data = reinterpret_cast<const double*>(face_data);
    normal[0] = *(normal_vector_data++);
    normal[1] = *(normal_vector_data++);
    normal[2] = *(normal_vector_data++);
    face_data = reinterpret_cast<const char*>(normal_vector_data);
    // M_surf matrix
    M_surf_data = reinterpret_cast<const double*>(face_data);
    face_data = reinterpret_cast<const char*>(M_surf_data + num_face_nodes * num_face_nodes);
    // IntS_shapeI_data
    IntS_shapeI_data = reinterpret_cast<const double*>(face_data);
    face_data = reinterpret_cast<const char*>(IntS_shapeI_data + num_face_nodes);
    // cell mapping data
    cell_mapping_data = reinterpret_cast<const std::uint32_t*>(face_data);
  }

  std::uint32_t num_face_nodes;
  double* outflow;
  std::array<double, 3> normal;
  const double* M_surf_data;
  const double* IntS_shapeI_data;
  const std::uint32_t* cell_mapping_data;
};

/// Cell view from contiguous block of memory
struct CellView
{
  __inline_host_dev__ CellView() {}
  __inline_host_dev__ void Update(const char* cell_data)
  {
    // number of faces and nodes (num_node = num_dof!)
    const std::uint32_t* num_node_and_face_data = reinterpret_cast<const std::uint32_t*>(cell_data);
    num_nodes = *(num_node_and_face_data++);
    num_faces = *(num_node_and_face_data++);
    cell_data = reinterpret_cast<const char*>(num_node_and_face_data);
    // density
    const double* density_data = reinterpret_cast<const double*>(cell_data);
    density = *(density_data++);
    cell_data = reinterpret_cast<const char*>(density_data);
    // total cross section pointer
    const double* const* total_xs_data = reinterpret_cast<const double* const*>(cell_data);
    total_xs = *(total_xs_data++);
    cell_data = reinterpret_cast<const char*>(total_xs_data);
    // phi address
    const std::uint64_t* phi_address_data = reinterpret_cast<const std::uint64_t*>(cell_data);
    phi_address = *(phi_address_data++);
    cell_data = reinterpret_cast<const char*>(phi_address_data);
    // save psi index
    const std::uint64_t* save_psi_index_data = reinterpret_cast<const std::uint64_t*>(cell_data);
    save_psi_index = *(save_psi_index_data++);
    cell_data = reinterpret_cast<const char*>(save_psi_index_data);
    // GM matrix
    GM_data = reinterpret_cast<const double*>(cell_data);
    cell_data = reinterpret_cast<const char*>(GM_data + num_nodes * num_nodes * 4);
    // face data
    offset_face_data = reinterpret_cast<const std::uint64_t*>(cell_data);
    face_data = reinterpret_cast<const char*>(offset_face_data + num_faces);
  }

  __inline_host_dev__ void GetFaceView(FaceView& face, const std::uint32_t& face_index)
  {
    face.Update(face_data + offset_face_data[face_index]);
  }

  std::uint32_t num_nodes;
  std::uint32_t num_faces;
  double density;
  const double* total_xs;
  std::uint64_t phi_address;
  std::uint64_t save_psi_index;
  const double* GM_data;
  const std::uint64_t* offset_face_data;
  const char* face_data;
};

/// Mesh view from contiguous block of memory
struct MeshView
{
  __inline_host_dev__ MeshView(const char* mesh_data)
  {
    const std::uint64_t* num_cells_data = reinterpret_cast<const std::uint64_t*>(mesh_data);
    num_cells = *(num_cells_data++);
    offset_cell_data = num_cells_data;
    cell_data = reinterpret_cast<const char*>(offset_cell_data + num_cells);
  }

  __inline_host_dev__ void GetCellView(CellView& cell, const std::uint32_t& cell_index)
  {
    cell.Update(cell_data + offset_cell_data[cell_index]);
  }

  const char* cell_data;
  std::uint64_t num_cells;
  const std::uint64_t* offset_cell_data;
};

} // namespace opensn
