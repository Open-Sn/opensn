// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <array>

namespace opensn
{
namespace fv_diffusion
{
class Boundary;

enum class BoundaryType : int
{
  Reflecting = 1,
  Dirichlet = 2,
  Neumann = 3,
  Robin = 4,
  Vacuum = 5
};

// ###################################################################
/**Parent class for diffusion boundaries*/
class Boundary
{
public:
  BoundaryType type_ = BoundaryType::Dirichlet;

  std::array<double, 3> values_ = {0.0, 0.0, 0.0};
};

} // namespace fv_diffusion
} // namespace opensn
