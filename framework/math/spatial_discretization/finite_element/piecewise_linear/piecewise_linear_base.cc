// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_base.h"
#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_slab_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_polygon_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_polyhedron_mapping.h"
#include "framework/mesh/mesh/mesh.h"

namespace opensn
{

PieceWiseLinearBase::PieceWiseLinearBase(const std::shared_ptr<Mesh>& grid,
                                         QuadratureOrder q_order,
                                         SpatialDiscretizationType sdm_type)
  : FiniteElementBase(grid, sdm_type, q_order),
    line_quad_order_arbitrary_(q_order),
    tri_quad_order_arbitrary_(q_order),
    tet_quad_order_arbitrary_(q_order)
{
}

void
PieceWiseLinearBase::CreateCellMappings()
{
  constexpr std::string_view fname = __PRETTY_FUNCTION__;

  auto MakeCellMapping = [&, this](auto cell_local_id) -> std::unique_ptr<CellMapping>
  {
    using namespace std;
    using namespace opensn;

    const auto& cell = grid_->GetLocalCell(cell_local_id);
    switch (cell.GetType())
    {
      case CellType::SLAB:
      {
        const auto& vol_quad = line_quad_order_arbitrary_;
        return make_unique<PieceWiseLinearSlabMapping>(cell_local_id, grid_, vol_quad);
      }
      case CellType::POLYGON:
      {
        const auto& vol_quad = tri_quad_order_arbitrary_;
        const auto& area_quad = line_quad_order_arbitrary_;
        return make_unique<PieceWiseLinearPolygonMapping>(
          cell_local_id, grid_, vol_quad, area_quad);
      }
      case CellType::POLYHEDRON:
      {
        const auto& vol_quad = tet_quad_order_arbitrary_;
        const auto& area_quad = tri_quad_order_arbitrary_;

        return make_unique<PieceWiseLinearPolyhedronMapping>(
          cell_local_id, grid_, vol_quad, area_quad);
      }
      default:
        throw std::invalid_argument(std::string(fname) +
                                    ": Unsupported cell type encountered. type_id=" +
                                    std::to_string(static_cast<int>(cell.GetType())));
    }
  };

  const auto ghost_ids = grid_->GetGhostGlobalIDs();
  cell_mappings_.reserve(grid_->GetLocalCellCount() + grid_->GhostCellCount());
  for (std::size_t cell_local_id = 0; cell_local_id < grid_->GetLocalCellCount(); ++cell_local_id)
    cell_mappings_.emplace_back(MakeCellMapping(cell_local_id));

  for (uint64_t ghost_id : ghost_ids)
  {
    auto cell_local_id = grid_->MapCellGlobalID2LocalID(ghost_id);
    auto ghost_mapping = MakeCellMapping(cell_local_id);
    cell_mappings_.emplace_back(std::move(ghost_mapping));
  }
}
} // namespace opensn
