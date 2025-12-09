// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

class SweepChunk;

enum class SchedulingAlgorithm
{
  FIRST_IN_FIRST_OUT = 1, ///< FIFO
  DEPTH_OF_GRAPH = 2      ///< DOG
};

class SweepScheduler
{
public:
  SweepScheduler(SchedulingAlgorithm scheduler_type,
                 AngleAggregation& angle_agg,
                 SweepChunk& sweep_chunk);

  void Sweep();

  void PrepareForSweep(bool use_boundary_source, bool zero_incoming_delayed_psi);

private:
  /// Applies a first-in-first-out sweep scheduling.
  void ScheduleAlgoFIFO(SweepChunk& sweep_chunk);

  /// Initializes the depth-of-graph algorithm.
  void InitializeAlgoDOG();

  /// Executes the depth-of-graph algorithm.
  void ScheduleAlgoDOG(SweepChunk& sweep_chunk);

private:
  SchedulingAlgorithm scheduler_type_;
  AngleAggregation& angle_agg_;
  SweepChunk& sweep_chunk_;

  struct RuleValues
  {
    std::shared_ptr<AngleSet> angle_set;
    int depth_of_graph;
    int sign_of_omegax;
    int sign_of_omegay;
    int sign_of_omegaz;
    size_t set_index;

    explicit RuleValues(std::shared_ptr<AngleSet>& ref_as)
      : angle_set(ref_as),
        depth_of_graph(0),
        sign_of_omegax(1),
        sign_of_omegay(1),
        sign_of_omegaz(1),
        set_index(0)
    {
    }
  };
  std::vector<RuleValues> rule_values_;
};

} // namespace opensn
