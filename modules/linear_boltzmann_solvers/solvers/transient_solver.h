// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include <functional>

namespace opensn
{

class DiscreteOrdinatesProblem;

class TransientKEigenSolver : public PowerIterationKEigenSolver
{
public:
  explicit TransientKEigenSolver(const InputParameters& params);

  ~TransientKEigenSolver() override = default;

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

private:
  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;

  std::vector<double>& phi_new_local_;
  std::vector<double>& precursor_new_local_;
  std::vector<std::vector<double>>& psi_new_local_;

  /// Previous time step vectors
  std::vector<double> phi_prev_local_;
  std::vector<double> precursor_prev_local_;
  std::vector<std::vector<double>> psi_prev_local_;

  /// Time discretization values and methods
  double stop_time_ = 0.1;
  double current_time_ = 0.0;
  unsigned int step_ = 0;
  bool verbose_ = true;
  bool initialized_ = false;
  std::function<void()> pre_advance_callback_;
  std::function<void()> post_advance_callback_;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<TransientKEigenSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
