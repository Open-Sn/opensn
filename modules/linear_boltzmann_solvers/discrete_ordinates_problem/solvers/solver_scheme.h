// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/source_functions/source_flags.h"
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace opensn
{

class AGSLinearSolver;
class DiscreteOrdinatesProblem;
class LBSGroupset;
class LinearSolver;
class SweepChunk;
struct WGSContext;

using SweepChunkFactory = std::function<std::shared_ptr<SweepChunk>(LBSGroupset&)>;

struct SolverScheme
{
  std::shared_ptr<AGSLinearSolver> ags_solver;
  std::vector<std::shared_ptr<WGSContext>> wgs_contexts;
  std::vector<std::shared_ptr<LinearSolver>> wgs_solvers;
  unsigned int max_groupset_size = 0;
  std::size_t max_level_size = 0;
  std::size_t max_angleset_size = 0;
};

SolverScheme BuildSolvers(DiscreteOrdinatesProblem& problem,
                          const SetSourceFunction& set_source_function,
                          const SweepChunkFactory& sweep_chunk_factory);

} // namespace opensn
