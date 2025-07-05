// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"

namespace opensn
{

MeshCarrier::MeshCarrier(LBSProblem& lbs_problem, TotalXSCarrier& xs, OutflowCarrier& outflow)
{
  std::uint64_t size = ComputeSize(lbs_problem);
  host_memory_.reserve(size);
  host_memory_.resize(size);
  Assemble(lbs_problem, xs, outflow);
  device_memory_ = crb::DeviceMemory<char>(size);
  crb::copy(device_memory_, host_memory_, size);
}

std::uint64_t
MeshCarrier::ComputeSize(LBSProblem& lbs_problem)
{
  std::uint64_t alloc_size = 0;
  // number of cells in the mesh
  alloc_size += sizeof(std::uint64_t);
  // offset of the data of each cell wrt to the origin pointer
  MeshContinuum& mesh = *(lbs_problem.GetGrid());
  const std::vector<UnitCellMatrices>& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  const SpatialDiscretization& discretization = lbs_problem.GetSpatialDiscretization();
  alloc_size += mesh.local_cells.size() * sizeof(std::uint64_t);
  // compute size for each cell in the mesh
  for (const Cell& cell : mesh.local_cells)
  {
    // number of faces and nodes
    alloc_size += 2 * sizeof(std::uint32_t);
    // densities
    alloc_size += sizeof(double);
    // pointer to total cross sections
    alloc_size += sizeof(std::uintptr_t);
    // phi address
    alloc_size += sizeof(std::uint64_t);
    // G and M matrix
    const UnitCellMatrices& unit_matrices = unit_cell_matrices[cell.local_id];
    const DenseMatrix<double>& M = unit_matrices.intV_shapeI_shapeJ;
    alloc_size += M.size() * (4 * sizeof(double));
    // offset to the data of each face
    std::size_t cell_num_faces = cell.faces.size();
    alloc_size += cell_num_faces * sizeof(std::uint64_t);
    // data of each face
    const std::vector<std::vector<int>>& face_node_mappings =
      discretization.GetCellMapping(cell).GetFaceNodeMappings();
    for (int f = 0; f < cell_num_faces; ++f)
    {
      // num_face_nodes
      alloc_size += sizeof(std::uint64_t);
      // pointer to outflow if applicable
      alloc_size += sizeof(std::uintptr_t);
      // normal vector
      alloc_size += 3 * sizeof(double);
      // M_surf matrix (num_face_nodes x num_face_nodes)
      std::size_t num_face_nodes = face_node_mappings[f].size();
      alloc_size += num_face_nodes * num_face_nodes * sizeof(double);
      // IntF_shape_I
      alloc_size += num_face_nodes * sizeof(double);
      // cell mapping (allocate up to the smallest even number to ensure 64-bit alignment)
      num_face_nodes = (num_face_nodes + 1) & ~static_cast<std::size_t>(1);
      alloc_size += num_face_nodes * sizeof(std::uint32_t);
    }
  }
  return alloc_size;
}

void
MeshCarrier::Assemble(LBSProblem& lbs_problem, TotalXSCarrier& xs, OutflowCarrier& outflow)
{
  // get information
  MeshContinuum& mesh = *(lbs_problem.GetGrid());
  const std::vector<UnitCellMatrices>& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  const SpatialDiscretization& discretization = lbs_problem.GetSpatialDiscretization();
  const std::vector<CellLBSView>& cell_transport_views = lbs_problem.GetCellTransportViews();
  // number of cells in the mesh
  char* data = host_memory_.data();
  std::uint64_t* num_cell_data = reinterpret_cast<std::uint64_t*>(data);
  std::uint64_t num_cells = mesh.local_cells.size();
  *(num_cell_data++) = num_cells;
  data = reinterpret_cast<char*>(num_cell_data);
  // copy data of each cell
  std::uint64_t* offset_cell_data = reinterpret_cast<std::uint64_t*>(data);
  data = reinterpret_cast<char*>(offset_cell_data + num_cells);
  for (char* cell_data = data; const Cell& cell : mesh.local_cells)
  {
    std::size_t cell_num_faces = cell.faces.size();
    const CellMapping& cell_mapping = discretization.GetCellMapping(cell);
    std::size_t cell_num_nodes = cell_mapping.GetNumNodes();
    // record current pointer offset
    *(offset_cell_data++) = cell_data - data;
    // number of faces and nodes (num_node = num_dof!)
    std::uint32_t* num_node_and_face_data = reinterpret_cast<std::uint32_t*>(cell_data);
    *(num_node_and_face_data++) = cell_num_nodes;
    *(num_node_and_face_data++) = cell_num_faces;
    cell_data = reinterpret_cast<char*>(num_node_and_face_data);
    // densitiy
    double* density_data = reinterpret_cast<double*>(cell_data);
    *(density_data++) = lbs_problem.GetDensitiesLocal()[cell.local_id];
    cell_data = reinterpret_cast<char*>(density_data);
    // pointer to total cross section
    double** total_xs_data = reinterpret_cast<double**>(cell_data);
    *(total_xs_data++) = xs.GetXSGPUData(cell.block_id);
    cell_data = reinterpret_cast<char*>(total_xs_data);
    // phi address
    std::uint64_t* phi_address_data = reinterpret_cast<std::uint64_t*>(cell_data);
    const CellLBSView& cell_transport_view = cell_transport_views[cell.local_id];
    std::uint64_t phi_address = cell_transport_view.MapDOF(0, 0, 0);
    *(phi_address_data++) = phi_address;
    cell_data = reinterpret_cast<char*>(phi_address_data);
    // GM matrices (G and M matrices are combined into one single Matrix of double4)
    const UnitCellMatrices& unit_matrices = unit_cell_matrices[cell.local_id];
    const DenseMatrix<Vector3>& G = unit_matrices.intV_shapeI_gradshapeJ;
    const DenseMatrix<double>& M = unit_matrices.intV_shapeI_shapeJ;
    double* GM_data = reinterpret_cast<double*>(cell_data);
    for (unsigned int i_row = 0; i_row < cell_num_nodes; ++i_row)
    {
      for (unsigned int i_col = 0; i_col < cell_num_nodes; ++i_col)
      {
        Vector3 element = G(i_row, i_col);
        *(GM_data++) = element.x;
        *(GM_data++) = element.y;
        *(GM_data++) = element.z;
        *(GM_data++) = M(i_row, i_col);
      }
    }
    cell_data = reinterpret_cast<char*>(GM_data);
    // face data
    std::uint64_t* offset_face_data = reinterpret_cast<std::uint64_t*>(cell_data);
    cell_data = reinterpret_cast<char*>(offset_face_data + cell_num_faces);
    const std::vector<DenseMatrix<double>>& M_surf_matrices = unit_matrices.intS_shapeI_shapeJ;
    const std::vector<Vector<double>>& IntS_shapeI_vectors = unit_matrices.intS_shapeI;
    const std::vector<std::vector<int>>& face_node_mappings = cell_mapping.GetFaceNodeMappings();
    char* face_data = cell_data;
    for (int f = 0; f < cell_num_faces; ++f)
    {
      const CellFace& face = cell.faces[f];
      *(offset_face_data++) = face_data - cell_data;
      // number of face node
      std::uint64_t* num_face_nodes_data = reinterpret_cast<std::uint64_t*>(face_data);
      std::size_t num_face_nodes = face_node_mappings[f].size();
      *(num_face_nodes_data++) = num_face_nodes;
      face_data = reinterpret_cast<char*>(num_face_nodes_data);
      // pointer to outflow
      double** outflow_data = reinterpret_cast<double**>(face_data);
      double* face_outflow = nullptr;
      bool is_boundary_face = not face.has_neighbor;
      if (is_boundary_face)
      {
        face_outflow = reinterpret_cast<double*>(outflow.GetDevicePtr());
        face_outflow += outflow.GetOffset(cell.local_id, f);
      }
      *(outflow_data++) = face_outflow;
      face_data = reinterpret_cast<char*>(outflow_data);
      // normal vector
      double* normal_vector_data = reinterpret_cast<double*>(face_data);
      *(normal_vector_data++) = face.normal.x;
      *(normal_vector_data++) = face.normal.y;
      *(normal_vector_data++) = face.normal.z;
      face_data = reinterpret_cast<char*>(normal_vector_data);
      // M_surf matrix (matrices are in row major)
      double* M_surf_data = reinterpret_cast<double*>(face_data);
      const DenseMatrix<double>& M_surf = M_surf_matrices[f];
      for (unsigned int fi = 0; fi < num_face_nodes; ++fi)
      {
        unsigned int i = cell_mapping.MapFaceNode(f, fi);
        for (unsigned int fj = 0; fj < num_face_nodes; ++fj)
        {
          unsigned int j = cell_mapping.MapFaceNode(f, fj);
          *(M_surf_data++) = M_surf(i, j);
        }
      }
      face_data = reinterpret_cast<char*>(M_surf_data);
      // IntS_shapeI
      double* IntS_shapeI_data = reinterpret_cast<double*>(face_data);
      const Vector<double>& IntS_shapeI = IntS_shapeI_vectors[f];
      for (unsigned int fi = 0; fi < num_face_nodes; ++fi)
      {
        unsigned int i = cell_mapping.MapFaceNode(f, fi);
        *(IntS_shapeI_data++) = IntS_shapeI(i);
      }
      face_data = reinterpret_cast<char*>(IntS_shapeI_data);
      // cell mapping data
      std::uint32_t* cell_mapping_data = reinterpret_cast<std::uint32_t*>(face_data);
      for (int i = 0; i < num_face_nodes; ++i)
      {
        cell_mapping_data[i] = face_node_mappings[f][i];
      }
      cell_mapping_data += ((num_face_nodes + 1) & ~1zu);
      face_data = reinterpret_cast<char*>(cell_mapping_data);
    }
    // put cell data pointer as face data end
    cell_data = face_data;
  }
}

} // namespace opensn
