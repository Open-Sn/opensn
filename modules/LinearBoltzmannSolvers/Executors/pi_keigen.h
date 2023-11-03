#pragma once

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"
#include "A_LBSSolver/IterativeMethods/wgs_context.h"

namespace lbs
{

class XXPowerIterationKEigen : public chi_physics::Solver
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
  std::shared_ptr<AGSLinearSolver<Mat, Vec, KSP>> primary_ags_solver_;
  lbs::SetSourceFunction active_set_source_function_;
  LBSGroupset& front_gs_;
  std::shared_ptr<chi_math::LinearSolver<Mat, Vec, KSP>> front_wgs_solver_;
  std::shared_ptr<lbs::WGSContext<Mat, Vec, KSP>> front_wgs_context_;

  double k_eff_ = 1.0;

public:
  static chi::InputParameters GetInputParameters();

  explicit XXPowerIterationKEigen(const chi::InputParameters& params);

  void Initialize() override;
  void Execute() override;

protected:
  void SetLBSFissionSource(const VecDbl& input, bool additive);
  void SetLBSScatterSource(const VecDbl& input, bool additive, bool suppress_wg_scat = false);
};

} // namespace lbs


