#include "opensn/framework/math/SpatialDiscretization/FiniteElement/FiniteElementBase.h"

chi_math::QuadratureOrder
chi_math::spatial_discretization::FiniteElementBase::GetQuadratureOrder() const
{
  return q_order_;
}
