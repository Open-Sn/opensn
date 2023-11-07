#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_base.h"

#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/math/spatial_discretization/cell_mappings/piecewise_linear/piecewise_linear_slab_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/piecewise_linear/piecewise_linear_polygon_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/piecewise_linear/piecewise_linear_polyhedron_mapping.h"

#define UnsupportedCellType(fname)                                                                 \
  std::invalid_argument((fname) + ": Unsupported cell type encountered. type_id=" +                \
                        std::to_string(static_cast<int>(cell.Type())));

namespace chi_math::spatial_discretization
{

PieceWiseLinearBase::PieceWiseLinearBase(const chi_mesh::MeshContinuum& grid,
                                         QuadratureOrder q_order,
                                         SDMType sdm_type,
                                         CoordinateSystemType cs_type)
  : FiniteElementBase(grid, cs_type, sdm_type, q_order),
    line_quad_order_arbitrary_(q_order),
    tri_quad_order_arbitrary_(q_order),
    quad_quad_order_arbitrary_(q_order),
    tet_quad_order_arbitrary_(q_order)
{
}

void
PieceWiseLinearBase::CreateCellMappings()
{
  constexpr std::string_view fname = __PRETTY_FUNCTION__;

  typedef cell_mapping::PieceWiseLinearSlabMapping SlabSlab;
  typedef cell_mapping::PieceWiseLinearPolygonMapping Polygon;
  typedef cell_mapping::PieceWiseLinearPolyhedronMapping Polyhedron;

  auto MakeCellMapping = [this, fname](const chi_mesh::Cell& cell)
  {
    using namespace std;
    using namespace chi_math;
    std::unique_ptr<chi_math::CellMapping> mapping;

    switch (cell.Type())
    {
      case chi_mesh::CellType::SLAB:
      {
        const auto& vol_quad = line_quad_order_arbitrary_;

        mapping = make_unique<SlabSlab>(cell, ref_grid_, vol_quad);
        break;
      }
      case chi_mesh::CellType::POLYGON:
      {
        const auto& vol_quad = tri_quad_order_arbitrary_;
        const auto& area_quad = line_quad_order_arbitrary_;

        mapping = make_unique<Polygon>(cell, ref_grid_, vol_quad, area_quad);
        break;
      }
      case chi_mesh::CellType::POLYHEDRON:
      {
        const auto& vol_quad = tet_quad_order_arbitrary_;
        const auto& area_quad = tri_quad_order_arbitrary_;

        mapping = make_unique<Polyhedron>(cell, ref_grid_, vol_quad, area_quad);
        break;
      }
      default:
        throw UnsupportedCellType(std::string(fname))
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
} // namespace chi_math::spatial_discretization
