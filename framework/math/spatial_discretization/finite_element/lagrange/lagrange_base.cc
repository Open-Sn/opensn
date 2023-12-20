#include "framework/math/spatial_discretization/finite_element/lagrange/lagrange_base.h"

#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_slab_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_quad_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_triangle_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_hex_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_wedge_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_tet_mapping.h"

#include "framework/math/spatial_discretization/cell_mappings/piecewise_linear/piecewise_linear_polygon_mapping.h"
#include "framework/math/spatial_discretization/cell_mappings/piecewise_linear/piecewise_linear_polyhedron_mapping.h"

#include "framework/logging/log_exceptions.h"

namespace opensn
{

LagrangeBase::LagrangeBase(const MeshContinuum& grid,
                           QuadratureOrder q_order,
                           SDMType sdm_type,
                           CoordinateSystemType cs_type)
  : FiniteElementBase(grid, cs_type, sdm_type, q_order),
    line_quad_order_arbitrary_(q_order),
    tri_quad_order_arbitrary_(q_order),
    quad_quad_order_arbitrary_(q_order),
    tet_quad_order_arbitrary_(q_order),
    hex_quad_order_arbitrary_(q_order),
    wedge_quad_order_arbitrary_(q_order)
{
  line_quad_order_arbitrary_.SetRange({-1.0, 1.0});

  CreateCellMappings();
}

void
LagrangeBase::CreateCellMappings()
{
  typedef LagrangeSlabMapping Slab;
  typedef LagrangeQuadMapping Quad;
  typedef LagrangeTriangleMapping Triangle;
  typedef LagrangeHexMapping Hex;
  typedef LagrangeWedgeMapping Wedge;
  typedef LagrangeTetMapping Tetrahedron;

  typedef PieceWiseLinearPolygonMapping Polygon;
  typedef PieceWiseLinearPolyhedronMapping Polyhedron;

  auto MakeCellMapping = [this](const Cell& cell)
  {
    using namespace std;
    using namespace opensn;
    std::unique_ptr<CellMapping> mapping;

    switch (cell.Type())
    {
      case CellType::SLAB:
      {
        const auto& vol_quad = line_quad_order_arbitrary_;
        const auto& area_quad = point_quadrature_;

        mapping = make_unique<Slab>(ref_grid_, cell, vol_quad, area_quad);
        break;
      }
      case CellType::POLYGON:
      {
        if (cell.SubType() == CellType::QUADRILATERAL)
        {
          const auto& vol_quad = quad_quad_order_arbitrary_;
          const auto& area_quad = line_quad_order_arbitrary_;

          mapping = make_unique<Quad>(ref_grid_, cell, vol_quad, area_quad);
          break;
        }
        else if (cell.SubType() == CellType::TRIANGLE)
        {
          const auto& vol_quad = tri_quad_order_arbitrary_;
          const auto& area_quad = line_quad_order_arbitrary_;

          mapping = make_unique<Triangle>(ref_grid_, cell, vol_quad, area_quad);
          break;
        }
        else
        {
          const auto& vol_quad = tri_quad_order_arbitrary_;
          const auto& area_quad = line_quad_order_arbitrary_;

          mapping = make_unique<Polygon>(cell, ref_grid_, vol_quad, area_quad);
          break;
        }
      }
      case CellType::POLYHEDRON:
      {
        if (cell.SubType() == CellType::HEXAHEDRON)
        {
          const auto& vol_quad = hex_quad_order_arbitrary_;
          const auto& area_quad = quad_quad_order_arbitrary_;

          mapping = make_unique<Hex>(ref_grid_, cell, vol_quad, area_quad);
          break;
        }
        else if (cell.SubType() == CellType::WEDGE)
        {
          const auto& vol_quad = wedge_quad_order_arbitrary_;
          const auto& area_quad1 = quad_quad_order_arbitrary_;
          const auto& area_quad2 = tri_quad_order_arbitrary_;

          mapping = make_unique<Wedge>(ref_grid_, cell, vol_quad, area_quad1, area_quad2);
          break;
        }
        else if (cell.SubType() == CellType::TETRAHEDRON)
        {
          const auto& vol_quad = tet_quad_order_arbitrary_;
          const auto& area_quad = tri_quad_order_arbitrary_;

          mapping = make_unique<Tetrahedron>(ref_grid_, cell, vol_quad, area_quad);
          break;
        }
        else
        {
          const auto& vol_quad = tet_quad_order_arbitrary_;
          const auto& area_quad = tri_quad_order_arbitrary_;

          mapping = make_unique<Polyhedron>(cell, ref_grid_, vol_quad, area_quad);
          break;
        }
      }
      default:
        ChiInvalidArgument("Unsupported cell type encountered");
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
