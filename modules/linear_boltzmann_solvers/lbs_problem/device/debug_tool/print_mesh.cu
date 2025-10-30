// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/device/debug_tool/print_mesh.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/mesh_view.h"
#include <cstddef>

#define BOX_DRAWING_VERTICAL "\xe2\x94\x82"
#define INDENT_0 "[DB] "
#define INDENT_1 "[DB]   "
#define INDENT_2 "[DB]     "
#define INDENT_3 "[DB]       "
#define INDENT_4 "[DB]         "
#define INDENT_5 "[DB]           "
#define INDENT(n) INDENT_##n

namespace opensn
{

namespace device_dbg
{

inline __device__ void
PrintFace(FaceView& face)
{
  // basic data
  std::printf(INDENT(4) "num_face_nodes: %" PRIu32 "\n", face.num_face_nodes);
  std::printf(INDENT(4) "outflow: %p %s\n",
              face.outflow,
              ((face.outflow) ? "-> boundary face" : "         -> non-boundary face"));
  std::printf(INDENT(4) "normal: (%f %f %f)\n", face.normal[0], face.normal[1], face.normal[2]);
  std::printf(INDENT(4) "M_surf:\n");
  // M_surf
  for (std::uint32_t i = 0; i < face.num_face_nodes; ++i)
  {
    std::printf(INDENT(5) "" BOX_DRAWING_VERTICAL);
    for (std::uint32_t j = 0; j < face.num_face_nodes; ++j)
    {
      std::printf("%s%8.6lf", ((j != 0) ? " " : ""), face.M_surf_data[i * face.num_face_nodes + j]);
    }
    std::printf(BOX_DRAWING_VERTICAL "\n");
  }
  // IntS_shapeI
  std::printf(INDENT(4) "IntS_shapeI: [");
  for (std::uint32_t i = 0; i < face.num_face_nodes; ++i)
  {
    std::printf("%s%8.6f", ((i != 0) ? " " : ""), face.IntS_shapeI_data[i]);
  }
  std::printf("]\n");
  // cell mapping
  std::printf(INDENT(4) "Cell mapping (face node -> cell node):\n");
  for (std::uint32_t i = 0; i < face.num_face_nodes; ++i)
  {
    std::printf(INDENT(5) "%" PRIu32 " -> %" PRIu32 "\n", i, face.cell_mapping_data[i]);
  }
}

inline __device__ void
PrintCell(CellView& cell)
{
  // basic data
  std::printf(INDENT(2) "num_nodes: %" PRIu32 "\n", cell.num_nodes);
  std::printf(INDENT(2) "num_faces: %" PRIu32 "\n", cell.num_faces);
  std::printf(INDENT(2) "density: %lf\n", cell.density);
  std::printf(INDENT(2) "total_xs_ptr: %p\n", cell.total_xs);
  std::printf(INDENT(2) "phi-address: %" PRIu64 "\n", cell.phi_address);
  std::printf(INDENT(2) "save_psi_index: %" PRIu64 "\n", cell.save_psi_index);
  // GM matrix
  std::printf(INDENT(2) "G Matrix:\n");
  for (std::uint32_t i = 0; i < cell.num_nodes; ++i)
  {
    std::printf(INDENT(3) "" BOX_DRAWING_VERTICAL);
    const double* GM_row = cell.GM_data + i * cell.num_nodes * 4;
    for (std::uint32_t j = 0; j < cell.num_nodes; ++j)
    {
      std::printf("%s(%8.6lf %8.6lf %8.6lf)",
                  ((j != 0) ? " " : ""),
                  GM_row[j * 4],
                  GM_row[j * 4 + 1],
                  GM_row[j * 4 + 2]);
    }
    std::printf(BOX_DRAWING_VERTICAL "\n");
  }
  std::printf(INDENT(2) "M Matrix:\n");
  for (std::uint32_t i = 0; i < cell.num_nodes; ++i)
  {
    std::printf(INDENT(3) BOX_DRAWING_VERTICAL);
    const double* GM_row = cell.GM_data + i * cell.num_nodes * 4;
    for (std::uint32_t j = 0; j < cell.num_nodes; ++j)
    {
      std::printf("%s%8.6lf", ((j != 0) ? " " : ""), GM_row[j * 4 + 3]);
    }
    std::printf(BOX_DRAWING_VERTICAL "\n");
  }
  // faces
  std::printf(INDENT(2) "Faces:\n");
  for (std::uint32_t f = 0; f < cell.num_faces; ++f)
  {
    std::printf(INDENT(3) "Face %" PRIu32 ":\n", f);
    FaceView face;
    cell.GetFaceView(face, f);
    PrintFace(face);
  }
}

__global__ void
PrintMesh(char* mesh_data)
{
  MeshView mesh(mesh_data);
  std::printf(INDENT(0) "Print mesh on GPU with %" PRIu64 " cells:\n", mesh.num_cells);
  for (std::uint32_t cell_idx = 0; cell_idx < mesh.num_cells; ++cell_idx)
  {
    std::printf(INDENT(1) "Cell %" PRIu32 ":\n", cell_idx);
    CellView cell;
    mesh.GetCellView(cell, cell_idx);
    PrintCell(cell);
  }
}

} // namespace device_dbg

void
device_dbg::PrintMeshOnDevice(MeshCarrier& mesh)
{
  device_dbg::PrintMesh<<<1, 1>>>(mesh.GetDevicePtr());
  crb::synchronize();
}

} // namespace opensn
