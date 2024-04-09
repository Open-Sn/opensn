// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/finite_element/finite_element_base.h"

namespace opensn
{

QuadratureOrder
FiniteElementBase::GetQuadratureOrder() const
{
  return q_order_;
}

} // namespace opensn
