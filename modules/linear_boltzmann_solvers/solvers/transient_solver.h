// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include <functional>

namespace opensn
{

class DiscreteOrdinatesProblem;
class SweepChunk;

class TransientKEigenSolver : public PowerIterationKEigenSolver
{
public:
  static InputParameters GetInputParameters();

  /// Options for initial condition normalization
  enum class NormalizationMethod
  {
    TOTAL_POWER = 0,    ///< Total reactor power
    POWER_DENSITY = 1,  ///< Power density
    NONE = 2            ///< No normalization
  };

  void SetPreAdvanceCallback(std::function<void()> callback);

  void SetPreAdvanceCallback(std::nullptr_t);

  void SetPostAdvanceCallback(std::function<void()> callback);

  void SetPostAdvanceCallback(std::nullptr_t);

  double GetCurrentTime() const { return current_time_; }

  unsigned int GetStep() const { return step_; }


protected:
  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;

  std::vector<double>& phi_new_local_;
  std::vector<double>& precursor_new_local_;
  std::vector<std::vector<double>>& psi_new_local_;

  /// Previous time step vectors
  std::vector<double> phi_prev_local_;
  std::vector<double> precursor_prev_local_;
  std::vector<std::vector<double>> psi_prev_local_;

  /// Fission rate vector
  std::vector<double> fission_rate_local_;

  /// Temporal domain and discretization information
  double dt_ = 2.0e-3;
  double theta_ = 0.5;
  double stop_time_ = 0.1;
  double current_time_ = 0.0;
  unsigned int step_ = 0;
  bool verbose_ = true;
  std::function<void()> pre_advance_callback_;
  std::function<void()> post_advance_callback_;

public:
  explicit TransientKEigenSolver(const InputParameters& params);

  static std::shared_ptr<TransientKEigenSolver> Create(const ParameterBlock& params);

  ~TransientKEigenSolver() override = default;

  void Initialize() override;
  void Execute() override;
  void Step();
  void Advance() override;
  void SetTimeStep(double dt);
  void SetTheta(double theta);

  //double ComputeBeta();
  void StepPrecursors();
};

}
