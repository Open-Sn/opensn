#include "framework/math/SpatialDiscretization/FiniteElement/Lagrange/LagrangeBase.h"

#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "framework/math/SpatialDiscretization/CellMappings/Lagrange/LagrangeSlabMapping.h"
#include "framework/math/SpatialDiscretization/CellMappings/Lagrange/LagrangeQuadMapping.h"
#include "framework/math/SpatialDiscretization/CellMappings/Lagrange/LagrangeTriangleMapping.h"
#include "framework/math/SpatialDiscretization/CellMappings/Lagrange/LagrangeHexMapping.h"
#include "framework/math/SpatialDiscretization/CellMappings/Lagrange/LagrangeWedgeMapping.h"
#include "framework/math/SpatialDiscretization/CellMappings/Lagrange/LagrangeTetMapping.h"

#include "framework/math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearPolygonMapping.h"
#include "framework/math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearPolyhedronMapping.h"

#include "framework/logging/chi_log_exceptions.h"

namespace chi_math::spatial_discretization
{

LagrangeBase::LagrangeBase(const chi_mesh::MeshContinuum& grid,
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
  typedef cell_mapping::LagrangeSlabMapping Slab;
  typedef cell_mapping::LagrangeQuadMapping Quad;
  typedef cell_mapping::LagrangeTriangleMapping Triangle;
  typedef cell_mapping::LagrangeHexMapping Hex;
  typedef cell_mapping::LagrangeWedgeMapping Wedge;
  typedef cell_mapping::LagrangeTetMapping Tetrahedron;

  typedef cell_mapping::PieceWiseLinearPolygonMapping Polygon;
  typedef cell_mapping::PieceWiseLinearPolyhedronMapping Polyhedron;

  auto MakeCellMapping = [this](const chi_mesh::Cell& cell)
  {
    using namespace std;
    using namespace chi_math;
    std::unique_ptr<chi_math::CellMapping> mapping;

    switch (cell.Type())
    {
      case chi_mesh::CellType::SLAB:
      {
        const auto& vol_quad = line_quad_order_arbitrary_;
        const auto& area_quad = point_quadrature_;

        mapping = make_unique<Slab>(ref_grid_, cell, vol_quad, area_quad);
        break;
      }
      case chi_mesh::CellType::POLYGON:
      {
        if (cell.SubType() == chi_mesh::CellType::QUADRILATERAL)
        {
          const auto& vol_quad = quad_quad_order_arbitrary_;
          const auto& area_quad = line_quad_order_arbitrary_;

          mapping = make_unique<Quad>(ref_grid_, cell, vol_quad, area_quad);
          break;
        }
        else if (cell.SubType() == chi_mesh::CellType::TRIANGLE)
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
      case chi_mesh::CellType::POLYHEDRON:
      {
        if (cell.SubType() == chi_mesh::CellType::HEXAHEDRON)
        {
          const auto& vol_quad = hex_quad_order_arbitrary_;
          const auto& area_quad = quad_quad_order_arbitrary_;

          mapping = make_unique<Hex>(ref_grid_, cell, vol_quad, area_quad);
          break;
        }
        else if (cell.SubType() == chi_mesh::CellType::WEDGE)
        {
          const auto& vol_quad = wedge_quad_order_arbitrary_;
          const auto& area_quad1 = quad_quad_order_arbitrary_;
          const auto& area_quad2 = tri_quad_order_arbitrary_;

          mapping = make_unique<Wedge>(ref_grid_, cell, vol_quad, area_quad1, area_quad2);
          break;
        }
        else if (cell.SubType() == chi_mesh::CellType::TETRAHEDRON)
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

} // namespace chi_math::spatial_discretization
