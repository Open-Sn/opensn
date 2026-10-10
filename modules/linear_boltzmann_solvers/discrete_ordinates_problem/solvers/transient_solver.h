// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"
#include <memory>
#include <string>
#include <functional>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

class TransientSolver : public DiscreteOrdinatesSolver
{
public:
  explicit TransientSolver(const InputParameters& params);

  ~TransientSolver() override = default;

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

protected:
  /// Requires time-dependent mode, except before Initialize() loads a steady-state initial
  /// condition (`read_initial_condition_path`), which switches the mode itself.
  void CheckRequirements(std::vector<std::string>& errors) const override;
  void InitializeSolver() override;
  void ExecuteSolver() override;
  void AdvanceSolver() override;

private:
  /// Sets this solver's time step and theta on the problem.
  void ApplyTimeParameters();

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

  /// Previous time step vectors
  std::vector<double> phi_prev_local_;
  std::vector<double> precursor_prev_local_;

  /// State tracked by CheckPrecursorStatus() to detect transitions.
  bool precursor_status_reported_ = false;
  bool last_use_precursors_ = false;
  bool last_has_fissionable_material_ = false;
  bool last_has_precursor_data_ = false;

  /// Time discretization values and methods. dt_ and theta_ are this solver's settings; the
  /// problem holds the values of the step being taken.
  double dt_ = 0.0;
  double theta_ = 0.0;
  double stop_time_ = 0.1;
  double current_time_ = 0.0;
  /// Time-discretization parameters used by the most recent completed Advance.
  double last_dt_ = 0.0;
  double last_theta_ = 0.0;
  unsigned int step_ = 0;
  bool verbose_ = true;
  bool enforce_stop_time_ = false;
  std::string initial_state_;
  std::function<void()> pre_advance_callback_;
  std::function<void()> post_advance_callback_;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<TransientSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
