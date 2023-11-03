#pragma once

#include "math/SpatialDiscretization/FiniteElement/FiniteElementBase.h"

#include "math/Quadratures/point_quadrature.h"
#include "math/Quadratures/quadrature_line.h"
#include "math/Quadratures/quadrature_triangle.h"
#include "math/Quadratures/quadrature_quadrilateral.h"
#include "math/Quadratures/quadrature_tetrahedron.h"
#include "math/Quadratures/quadrature_hexahedron.h"
#include "math/Quadratures/quadrature_wedge.h"

namespace chi_math::spatial_discretization
{

/**Base class for Lagrange spatial discretizations.
 * \ingroup doc_SpatialDiscretization*/
class LagrangeBase : public FiniteElementBase
{
protected:
  LagrangeBase(const chi_mesh::MeshContinuum& grid,
               QuadratureOrder q_order,
               SDMType sdm_type,
               CoordinateSystemType cs_type);

  PointQuadrature point_quadrature_;
  QuadratureLine line_quad_order_arbitrary_;
  QuadratureTriangle tri_quad_order_arbitrary_;
  QuadratureQuadrilateral quad_quad_order_arbitrary_;
  QuadratureTetrahedron tet_quad_order_arbitrary_;
  QuadratureHexahedron hex_quad_order_arbitrary_;
  QuadratureWedge wedge_quad_order_arbitrary_;

  void CreateCellMappings();
};

} // namespace chi_math::spatial_discretization
