// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/finite_element/finite_element_base.h"
#include "framework/math/quadratures/spatial/line_quadrature.h"
#include "framework/math/quadratures/spatial/triangle_quadrature.h"
#include "framework/math/quadratures/spatial/tetrahedra_quadrature.h"

namespace opensn
{

/**
 * Base class for PieceWiseLinear based discretization.
 *
 * \ingroup doc_SpatialDiscretization
 */
class PieceWiseLinearBase : public FiniteElementBase
{
protected:
  /// Constructor
  explicit PieceWiseLinearBase(const MeshContinuum& grid,
                               QuadratureOrder q_order,
                               SpatialDiscretizationType sdm_type,
                               CoordinateSystemType cs_type);

  LineQuadrature line_quad_order_arbitrary_;
  TriangleQuadrature tri_quad_order_arbitrary_;
  TetrahedraQuadrature tet_quad_order_arbitrary_;

  void CreateCellMappings();
};

} // namespace opensn
