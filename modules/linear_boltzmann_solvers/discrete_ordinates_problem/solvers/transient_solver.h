// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"
#include "framework/object_factory.h"
#include <memory>
#include <string>
#include <functional>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

class TransientSolver : public Solver
{
public:
  explicit TransientSolver(const InputParameters& params);

  ~TransientSolver() override = default;

  void Initialize() override;
  void Execute() override;
  void Advance() override;
  void SetTimeStep(double dt);
  void SetTheta(double theta);
  void StepPrecursors();
  void SetPreAdvanceCallback(std::function<void()> callback);
  void SetPreAdvanceCallback(std::nullptr_t);
  void SetPostAdvanceCallback(std::function<void()> callback);
  void SetPostAdvanceCallback(std::nullptr_t);
  double GetCurrentTime() const { return current_time_; }
  unsigned int GetStep() const { return step_; }
  BalanceTable ComputeBalanceTable() const;

private:
  bool ReadRestartData();
  bool ReadInitialConditionData();
  bool WriteRestartData();

  /**
   * Re-evaluates the current precursor state (use_precursors vs. presence of fissionable
   * material and delayed-neutron data in the active cross-section map) and logs a warning
   * or informational message only when that status has changed since the last check. This
   * keeps the reported status accurate across cross-section swaps (e.g. via
   * Problem.SetXSMap), rather than reflecting only a one-time check made at Initialize().
   */
  void CheckPrecursorStatus();

  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;

  /// Previous time step vectors
  std::vector<double> phi_prev_local_;
  std::vector<double> precursor_prev_local_;

  /// State tracked by CheckPrecursorStatus() to detect transitions.
  bool precursor_status_reported_ = false;
  bool last_use_precursors_ = false;
  bool last_has_fissionable_material_ = false;
  bool last_has_precursor_data_ = false;

  /// Time discretization values and methods
  double stop_time_ = 0.1;
  double current_time_ = 0.0;
  unsigned int step_ = 0;
  bool verbose_ = true;
  bool initialized_ = false;
  bool enforce_stop_time_ = false;
  std::string initial_state_;
  std::function<void()> pre_advance_callback_;
  std::function<void()> post_advance_callback_;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<TransientSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
