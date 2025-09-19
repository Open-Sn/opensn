// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"
#include "framework/math/math_time_stepping.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"

namespace opensn
{

class DiscreteOrdinatesProblem;
class SweepChunk;

class TransientKEigenSolver : public PowerIterationKEigenSolver
{
public:
  static InputParameters GetInputParameters();
  SteppingMethod method = SteppingMethod::CRANK_NICOLSON;

  /// Options for initial condition normalization
  enum class NormalizationMethod
  {
    TOTAL_POWER = 0,    ///< Total reactor power
    POWER_DENSITY = 1,  ///< Power density
    NONE = 2            ///< No normalization
  };

  struct Options
  {
    int verbosity_level = 1;

    bool inhibit_advance = false;
    double t_final = 0.1;
    int max_time_steps = 10;

    bool scale_fission_xs = false;
    NormalizationMethod normalization_method = NormalizationMethod::TOTAL_POWER;
  } transient_options_;


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
  double time_ = 0.0;

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

  //Iterative operations
  std::shared_ptr<SweepChunk> SetTransientSweepChunk(LBSGroupset& groupset);

  //double ComputeBeta();
  void PostStepCallBackFunction() const;
  void StepPrecursors();
};

}
