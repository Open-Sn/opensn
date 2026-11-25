// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"

namespace opensn
{

class DiscreteOrdinatesProblem;
class DiscreteOrdinatesKEigenAcceleration;

class PowerIterationKEigenSolver : public DiscreteOrdinatesSolver
{
public:
  explicit PowerIterationKEigenSolver(const InputParameters& params);
  /// Return the current k-eigenvalue.
  double GetEigenvalue() const { return k_eff_; }
  /// Return the k-eigenvalue convergence tolerance.
  double GetKTolerance() const { return k_tolerance_; }
  /// Return the number of completed power iterations.
  unsigned int GetNumPowerIterations() const { return num_power_iterations_; }
  /// Return true if the last Execute() met its convergence criteria before max_iters.
  bool IsConverged() const { return converged_; }
  /// Return the fission-source tolerance. Zero means the fission-source test is disabled.
  double GetFissionSourceTolerance() const { return fission_source_tolerance_; }
  /**
   * Return the root-mean-square relative change of the normalized cell fission sources in the
   * last power iteration. Zero when the fission-source test is disabled.
   */
  double GetFissionSourceChange() const { return fission_source_change_; }
  /**
   * Return the number of transport sweeps between acceleration corrections: sweeps_per_update
   * when it is set, otherwise 1.
   */
  unsigned int GetSweepsPerUpdate() const
  {
    return sweeps_per_update_ > 0 ? sweeps_per_update_ : 1;
  }
  /// Return true if sweeps_per_update was left to the automatic choice.
  bool IsSweepsPerUpdateAutomatic() const { return sweeps_per_update_ == 0; }

  BalanceTable ComputeBalanceTable() const;

  /// Combines function calls to set fission source.
  void SetLBSFissionSource(const std::vector<double>& input, bool additive);

  /// Combines function calls to set scattering source source.
  void SetLBSScatterSource(const std::vector<double>& input,
                           bool additive,
                           bool suppress_wg_scat = false);

protected:
  void CheckRequirements(std::vector<std::string>& errors) const override;
  void InitializeSolver() override;
  void ExecuteSolver() override;

  const std::shared_ptr<DiscreteOrdinatesKEigenAcceleration> acceleration_;

  unsigned int max_iters_;
  double k_eff_;
  double k_tolerance_;
  double fission_source_tolerance_;
  double fission_source_change_ = 0.0;
  /// Transport sweeps between acceleration corrections. Zero selects the automatic choice.
  const unsigned int sweeps_per_update_;
  double F_prev_;
  bool reset_phi0_;
  unsigned int num_power_iterations_ = 0;
  bool converged_ = false;

  std::vector<double>& q_moments_local_;
  std::vector<double>& phi_old_local_;
  std::vector<double>& phi_new_local_;

  const std::vector<LBSGroupset>& groupsets_;
  SetSourceFunction active_set_source_function_;

private:
  /// Refreshes the problem objects the solver uses and removes fission from the WGS operators
  /// until the returned guard is destroyed.
  [[nodiscard]] ScopedRestore PrepareProblem();

  bool WriteRestartData();

  bool ReadRestartData();

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<PowerIterationKEigenSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
