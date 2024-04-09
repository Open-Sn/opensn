// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <vector>

namespace opensn
{
namespace mg_diffusion
{
class Boundary;

enum class BoundaryType : int
{
  Reflecting = 1,
  Neumann = 3,
  Robin = 4,
  Vacuum = 5
};

// ###################################################################
/**Parent class for multigroup diffusion boundaries*/
class Boundary
{
public:
  BoundaryType type_ = BoundaryType::Vacuum;

  std::array<std::vector<double>, 3> mg_values_;
  // std::array<double, 3> mg_values = {0.25,0.5,0.};
};

} // namespace mg_diffusion
} // namespace opensn
