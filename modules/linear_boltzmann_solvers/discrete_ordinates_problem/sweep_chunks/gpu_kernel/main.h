// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/buffer.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include <utility>

namespace opensn::gpu_kernel
{

template <std::uint32_t... I, class F>
__device__ inline void
ForDOFs1ToN(std::integer_sequence<std::uint32_t, I...>, F&& f)
{
  (f(std::integral_constant<std::uint32_t, I + 1>{}), ...);
}

template <class F>
__device__ inline void
ForDOFs1ToMax(F&& f)
{
  ForDOFs1ToN(std::make_integer_sequence<std::uint32_t, LBSProblem::max_dofs_gpu>{},
              std::forward<F>(f));
}

template <class... Args>
__device__ inline void
SweepDispatch(std::uint32_t n, Args&&... args)
{
  bool done = false;
  ForDOFs1ToMax(
    [&](auto dof_c)
    {
      constexpr std::uint32_t dof = decltype(dof_c)::value;
      if (!done && n == dof)
      {
        gpu_kernel::Sweep<dof>(std::forward<Args>(args)...);
        done = true;
      }
    });
}

} // namespace opensn::gpu_kernel
