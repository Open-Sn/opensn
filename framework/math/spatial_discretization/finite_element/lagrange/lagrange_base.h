#pragma once

#include "framework/math/spatial_discretization/finite_element/finite_element_base.h"

#include "framework/math/quadratures/point_quadrature.h"
#include "framework/math/quadratures/quadrature_line.h"
#include "framework/math/quadratures/quadrature_triangle.h"
#include "framework/math/quadratures/quadrature_quadrilateral.h"
#include "framework/math/quadratures/quadrature_tetrahedron.h"
#include "framework/math/quadratures/quadrature_hexahedron.h"
#include "framework/math/quadratures/quadrature_wedge.h"

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
