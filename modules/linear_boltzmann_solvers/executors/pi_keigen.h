// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"

namespace opensn
{
namespace lbs
{

class PowerIterationKEigen : public opensn::Solver
{
protected:
  LBSSolver& lbs_solver_;

  size_t max_iters_;
  double k_eff_;
  double k_tolerance_;
  double F_prev_;
  bool reset_phi0_;

  std::vector<double>& q_moments_local_;
  std::vector<double>& phi_old_local_;
  std::vector<double>& phi_new_local_;

  std::vector<LBSGroupset>& groupsets_;
  std::shared_ptr<AGSLinearSolver> primary_ags_solver_;
  lbs::SetSourceFunction active_set_source_function_;

  LBSGroupset& front_gs_;
  std::shared_ptr<LinearSolver> front_wgs_solver_;
  std::shared_ptr<lbs::WGSContext> front_wgs_context_;

public:
  static InputParameters GetInputParameters();

  explicit PowerIterationKEigen(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

protected:
  /**
   * Combines function calls to set fission source.
   */
  void SetLBSFissionSource(const std::vector<double>& input, bool additive);

  /**
   * Combines function calls to set scattering source source.
   */
  void SetLBSScatterSource(const std::vector<double>& input,
                           bool additive,
                           bool suppress_wg_scat = false);

  void WriteRestartData();

  void ReadRestartData();
};

} // namespace lbs
} // namespace opensn
