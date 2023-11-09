#pragma once

#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/spatial_discretization/finite_element/quadrature_point_data.h"

namespace chi_math::spatial_discretization
{
/**Base Finite Element spatial discretization class.
 * \ingroup doc_SpatialDiscretization*/
class FiniteElementBase : public chi_math::SpatialDiscretization
{
public:
  QuadratureOrder GetQuadratureOrder() const;

protected:
  explicit FiniteElementBase(const chi_mesh::MeshContinuum& grid,
                             CoordinateSystemType cs_type,
                             SDMType sdm_type,
                             QuadratureOrder q_order)
    : SpatialDiscretization(grid, cs_type, sdm_type), q_order_(q_order)
  {
  }

  const QuadratureOrder q_order_;
};
} // namespace chi_math::spatial_discretization
