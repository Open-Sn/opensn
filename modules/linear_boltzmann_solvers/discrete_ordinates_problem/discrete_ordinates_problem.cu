// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "caribou/caribou.h"
namespace crb = caribou;

namespace opensn
{

void
DiscreteOrdinatesProblem::CheckSystem()
{
  std::uint32_t num_gpus = crb::get_num_gpus();
  if (num_gpus == 0)
  {
    throw std::runtime_error("No GPU detected.\n");
  }
}

} // namespace opensn
