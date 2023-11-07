#include "framework/math/spatial_discretization/finite_element/finite_element_base.h"

chi_math::QuadratureOrder
chi_math::spatial_discretization::FiniteElementBase::GetQuadratureOrder() const
{
  return q_order_;
}
