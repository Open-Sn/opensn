#pragma once

#include "framework/physics/solver_base/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/nl_keigen_ags_solver.h"

#include <petscsnes.h>

namespace opensn
{
namespace lbs
{

class XXNonLinearKEigen : public opensn::Solver
{
protected:
  LBSSolver& lbs_solver_;
  std::shared_ptr<NLKEigenAGSContext> nl_context_;
  NLKEigenvalueAGSSolver nl_solver_;

  bool reinit_phi_1_;
  int num_free_power_its_;

public:
  static InputParameters GetInputParameters();
  explicit XXNonLinearKEigen(const InputParameters& params);

  void Initialize() override;
  void Execute() override;
};

} // namespace lbs
} // namespace opensn
