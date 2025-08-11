// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"

namespace opensn
{

class DiscreteOrdinatesProblem;
class DiscreteOrdinatesKEigenAcceleration;
class AGSSolver;
class LinearSolver;

class PowerIterationKEigenSolver : public Solver
{
protected:
  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;
  const std::shared_ptr<DiscreteOrdinatesKEigenAcceleration> acceleration_;

  size_t max_iters_;
  double k_eff_;
  double k_tolerance_;
  double F_prev_;
  bool reset_phi0_;

  std::vector<double>& q_moments_local_;
  std::vector<double>& phi_old_local_;
  std::vector<double>& phi_new_local_;

  std::vector<LBSGroupset>& groupsets_;
  std::shared_ptr<AGSSolver> ags_solver_;
  SetSourceFunction active_set_source_function_;

  LBSGroupset& front_gs_;
  std::shared_ptr<LinearSolver> front_wgs_solver_;
  std::shared_ptr<WGSContext> front_wgs_context_;

public:
  explicit PowerIterationKEigenSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;
  /// Return the current k-eigenvalue
  double GetEigenvalue() const { return k_eff_; }

  /// Combines function calls to set fission source.
  void SetLBSFissionSource(const std::vector<double>& input, bool additive);

  /// Combines function calls to set scattering source source.
  void SetLBSScatterSource(const std::vector<double>& input,
                           bool additive,
                           bool suppress_wg_scat = false);

private:
  bool WriteRestartData();

  bool ReadRestartData();

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<PowerIterationKEigenSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
