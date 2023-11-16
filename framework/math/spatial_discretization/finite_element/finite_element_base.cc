#include "framework/math/spatial_discretization/finite_element/finite_element_base.h"

namespace opensn
{

QuadratureOrder
FiniteElementBase::GetQuadratureOrder() const
{
  return q_order_;
}

} // namespace opensn
