#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_base.h"
#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_slab_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_polygon_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_polyhedron_mapping.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

PieceWiseLinearBase::PieceWiseLinearBase(const MeshContinuum& grid,
                                         QuadratureOrder q_order,
                                         SDMType sdm_type,
                                         CoordinateSystemType cs_type)
  : FiniteElementBase(grid, cs_type, sdm_type, q_order),
    line_quad_order_arbitrary_(q_order),
    tri_quad_order_arbitrary_(q_order),
    tet_quad_order_arbitrary_(q_order)
{
}

void
PieceWiseLinearBase::CreateCellMappings()
{
  constexpr std::string_view fname = __PRETTY_FUNCTION__;

  typedef PieceWiseLinearSlabMapping SlabSlab;
  typedef PieceWiseLinearPolygonMapping Polygon;
  typedef PieceWiseLinearPolyhedronMapping Polyhedron;

  auto MakeCellMapping = [this, fname](const Cell& cell)
  {
    using namespace std;
    using namespace opensn;
    std::unique_ptr<CellMapping> mapping;

    switch (cell.Type())
    {
      case CellType::SLAB:
      {
        const auto& vol_quad = line_quad_order_arbitrary_;

        mapping = make_unique<SlabSlab>(cell, ref_grid_, vol_quad);
        break;
      }
      case CellType::POLYGON:
      {
        const auto& vol_quad = tri_quad_order_arbitrary_;
        const auto& area_quad = line_quad_order_arbitrary_;

        mapping = make_unique<Polygon>(cell, ref_grid_, vol_quad, area_quad);
        break;
      }
      case CellType::POLYHEDRON:
      {
        const auto& vol_quad = tet_quad_order_arbitrary_;
        const auto& area_quad = tri_quad_order_arbitrary_;

        mapping = make_unique<Polyhedron>(cell, ref_grid_, vol_quad, area_quad);
        break;
      }
      default:
        throw std::invalid_argument(std::string(fname) +
                                    ": Unsupported cell type encountered. type_id=" +
                                    std::to_string(static_cast<int>(cell.Type())));
    }
    return mapping;
  };

  for (const auto& cell : ref_grid_.local_cells)
    cell_mappings_.push_back(MakeCellMapping(cell));

  const auto ghost_ids = ref_grid_.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto ghost_mapping = MakeCellMapping(ref_grid_.cells[ghost_id]);
    nb_cell_mappings_.insert(std::make_pair(ghost_id, std::move(ghost_mapping)));
  }
}
} // namespace opensn
