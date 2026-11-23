// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_cell_sweep_data.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh/mesh.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_view.h"

namespace opensn
{

AAHCellSweepData::AAHCellSweepData(const Mesh& grid,
                                   const SpatialDiscretization& discretization,
                                   const std::vector<UnitCellMatrices>& unit_cell_matrices,
                                   const std::vector<CellLBSView>& cell_transport_views)
{
  const auto num_cells = grid.GetLocalCellCount();

  // Size the arrays first, so that they are allocated once.
  std::size_t num_faces = 0, num_vectors = 0, num_scalars = 0, num_face_nodes = 0;
  for (std::size_t c = 0; c < num_cells; ++c)
  {
    const auto& cell = grid.GetLocalCell(c);
    const auto& mapping = discretization.GetCellMapping(cell);
    const std::size_t n = mapping.GetNumNodes();
    num_faces += cell.faces.size();
    num_vectors += n * n;
    num_scalars += n * n;
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const std::size_t nf = mapping.GetNumFaceNodes(f);
      num_scalars += nf * nf + nf;
      num_face_nodes += nf;
    }
  }
  cells_.reserve(num_cells);
  faces_.reserve(num_faces);
  vectors_.reserve(num_vectors);
  scalars_.reserve(num_scalars);
  face_nodes_.reserve(num_face_nodes);

  for (std::size_t c = 0; c < num_cells; ++c)
  {
    const auto& cell = grid.GetLocalCell(c);
    const auto& mapping = discretization.GetCellMapping(cell);
    const auto& matrices = unit_cell_matrices[c];
    const auto& transport_view = cell_transport_views[c];
    const std::size_t n = mapping.GetNumNodes();

    Cell entry;
    entry.first_node = transport_view.GetFirstNode();
    entry.face_offset = faces_.size();
    entry.vector_offset = vectors_.size();
    entry.scalar_offset = scalars_.size();
    entry.num_nodes = static_cast<std::uint32_t>(n);
    entry.num_faces = static_cast<std::uint32_t>(cell.faces.size());
    cells_.push_back(entry);

    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        vectors_.push_back(matrices.intV_shapeI_gradshapeJ(i, j));
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        scalars_.push_back(matrices.intV_shapeI_shapeJ(i, j));

    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& cell_face = cell.faces[f];
      const std::size_t nf = mapping.GetNumFaceNodes(f);

      Face face;
      face.normal = cell_face.normal;
      face.neighbor_id = cell_face.neighbor_id;
      face.scalar_offset = scalars_.size();
      face.node_offset = face_nodes_.size();
      face.num_nodes = static_cast<std::uint32_t>(nf);
      face.has_neighbor = cell_face.has_neighbor;
      face.is_local = transport_view.IsFaceLocal(f);
      faces_.push_back(face);

      for (std::size_t fi = 0; fi < nf; ++fi)
        face_nodes_.push_back(mapping.MapFaceNode(f, fi));
      const auto& face_mass = matrices.intS_shapeI_shapeJ[f];
      for (std::size_t fi = 0; fi < nf; ++fi)
        for (std::size_t fj = 0; fj < nf; ++fj)
          scalars_.push_back(face_mass(mapping.MapFaceNode(f, fi), mapping.MapFaceNode(f, fj)));
      const auto& face_shape = matrices.intS_shapeI[f];
      for (std::size_t fi = 0; fi < nf; ++fi)
        scalars_.push_back(face_shape(mapping.MapFaceNode(f, fi)));
    }
  }
}

} // namespace opensn
