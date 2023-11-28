#pragma once

#include "framework/physics/solver_base/solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/nl_keigen_ags_solver.h"

#include <petscsnes.h>

namespace lbs
{

class XXNonLinearKEigen : public chi_physics::Solver
{
protected:
  LBSSolver& lbs_solver_;
  std::shared_ptr<NLKEigenAGSContext> nl_context_;
  NLKEigenvalueAGSSolver nl_solver_;

  bool reinit_phi_1_;
  int num_free_power_its_;

public:
  static chi::InputParameters GetInputParameters();
  explicit XXNonLinearKEigen(const chi::InputParameters& params);

  void Initialize() override;
  void Execute() override;
};

} // namespace lbs
