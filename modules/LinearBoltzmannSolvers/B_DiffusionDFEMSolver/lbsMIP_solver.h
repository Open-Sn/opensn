#pragma once

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

namespace lbs
{

class DiffusionDFEMSolver : public LBSSolver
{
protected:
  typedef std::shared_ptr<acceleration::DiffusionMIPSolver> MIPSolverPtr;

public:
  std::vector<MIPSolverPtr> gs_mip_solvers_;

public:
  // 00
  static chi::InputParameters GetInputParameters();
  explicit DiffusionDFEMSolver(const chi::InputParameters& params);
  ~DiffusionDFEMSolver() override;
  void Initialize() override;
  void InitializeWGSSolvers() override;
};

} // namespace lbs
