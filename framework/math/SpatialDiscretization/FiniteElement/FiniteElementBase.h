#ifndef SPATIAL_DISCRETIZATION_FE_H
#define SPATIAL_DISCRETIZATION_FE_H

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/UnknownManager/unknown_manager.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

// ###################################################################
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

#endif