// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include "framework/object_factory.h"
#include <memory>
#include <string>
#include <functional>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;
class AGSLinearSolver;

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

private:
  void RefreshLocalViews();
  void UpdateHasFissionableMaterial();

  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;
  std::shared_ptr<AGSLinearSolver> ags_solver_;

  std::vector<double>* q_moments_local_ = nullptr;
  std::vector<double>* phi_old_local_ = nullptr;

  std::vector<double>* phi_new_local_ = nullptr;
  std::vector<double>* precursor_new_local_ = nullptr;
  std::vector<std::vector<double>>* psi_new_local_ = nullptr;

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
  bool enforce_stop_time_ = false;
  bool has_fissionable_material_ = false;
  std::string initial_state_;
  std::function<void()> pre_advance_callback_;
  std::function<void()> post_advance_callback_;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<TransientSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
