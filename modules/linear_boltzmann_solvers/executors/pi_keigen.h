// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"

namespace opensn
{
namespace lbs
{

class XXPowerIterationKEigen : public opensn::Solver
{
protected:
  LBSSolver& lbs_solver_;
  size_t max_iters_;
  double k_tolerance_;
  bool reinit_phi_1_;

  VecDbl& q_moments_local_;
  VecDbl& phi_old_local_;
  VecDbl& phi_new_local_;
  std::vector<LBSGroupset>& groupsets_;
  std::shared_ptr<AGSLinearSolver> primary_ags_solver_;
  lbs::SetSourceFunction active_set_source_function_;
  LBSGroupset& front_gs_;
  std::shared_ptr<LinearSolver> front_wgs_solver_;
  std::shared_ptr<lbs::WGSContext> front_wgs_context_;

  double k_eff_ = 1.0;

public:
  static InputParameters GetInputParameters();

  explicit XXPowerIterationKEigen(const InputParameters& params);

  void Initialize() override;
  void Execute() override;

protected:
  /**
   * Combines function calls to set fission source.
   */
  void SetLBSFissionSource(const VecDbl& input, bool additive);

  /**
   * Combines function calls to set scattering source source.
   */
  void SetLBSScatterSource(const VecDbl& input, bool additive, bool suppress_wg_scat = false);
};

} // namespace lbs
} // namespace opensn
