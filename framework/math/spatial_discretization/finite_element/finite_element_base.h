// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/math/unknown_manager/unknown_manager.h"

namespace opensn
{

/**
 * Base Finite Element spatial discretization class.
 *
 * \ingroup doc_SpatialDiscretization
 */
class FiniteElementBase : public SpatialDiscretization
{
public:
  QuadratureOrder GetQuadratureOrder() const;

protected:
  explicit FiniteElementBase(const std::shared_ptr<MeshContinuum> grid,
                             SpatialDiscretizationType sdm_type,
                             QuadratureOrder q_order)
    : SpatialDiscretization(grid, sdm_type), q_order_(q_order)
  {
  }

  const QuadratureOrder q_order_;
};

} // namespace opensn
