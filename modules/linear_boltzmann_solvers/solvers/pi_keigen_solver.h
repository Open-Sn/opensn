// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"

namespace opensn
{

class DiscreteOrdinatesProblem;
class DiscreteOrdinatesKEigenAcceleration;

class PowerIterationKEigenSolver : public Solver
{
public:
  explicit PowerIterationKEigenSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;
  /// Return the current k-eigenvalue
  double GetEigenvalue() const { return k_eff_; }

  BalanceTable ComputeBalanceTable() const;

  /// Combines function calls to set fission source.
  void SetLBSFissionSource(const std::vector<double>& input, bool additive);

  /// Combines function calls to set scattering source source.
  void SetLBSScatterSource(const std::vector<double>& input,
                           bool additive,
                           bool suppress_wg_scat = false);

protected:
  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;
  const std::shared_ptr<DiscreteOrdinatesKEigenAcceleration> acceleration_;

  unsigned int max_iters_;
  double k_eff_;
  double k_tolerance_;
  double F_prev_;
  bool reset_phi0_;

  std::vector<double>& q_moments_local_;
  std::vector<double>& phi_old_local_;
  std::vector<double>& phi_new_local_;

  const std::vector<LBSGroupset>& groupsets_;
  SetSourceFunction active_set_source_function_;

  bool initialized_ = false;

private:
  bool WriteRestartData();

  bool ReadRestartData();

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<PowerIterationKEigenSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
