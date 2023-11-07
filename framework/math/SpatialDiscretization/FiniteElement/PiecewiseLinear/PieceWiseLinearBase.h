#pragma once

#include "framework/math/SpatialDiscretization/FiniteElement/FiniteElementBase.h"

#include "framework/math/Quadratures/quadrature_line.h"
#include "framework/math/Quadratures/quadrature_triangle.h"
#include "framework/math/Quadratures/quadrature_quadrilateral.h"
#include "framework/math/Quadratures/quadrature_tetrahedron.h"

namespace chi_math::spatial_discretization
{

/**Base class for PieceWiseLinear based discretization.
 * \ingroup doc_SpatialDiscretization*/
class PieceWiseLinearBase : public FiniteElementBase
{
protected:
  /**Constructor*/
  explicit PieceWiseLinearBase(const chi_mesh::MeshContinuum& grid,
                               QuadratureOrder q_order,
                               SDMType sdm_type,
                               CoordinateSystemType cs_type);

  QuadratureLine line_quad_order_arbitrary_;
  QuadratureTriangle tri_quad_order_arbitrary_;
  QuadratureQuadrilateral quad_quad_order_arbitrary_;
  QuadratureTetrahedron tet_quad_order_arbitrary_;

  void CreateCellMappings();
};

} // namespace chi_math::spatial_discretization
