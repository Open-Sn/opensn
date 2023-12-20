#pragma once

#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

namespace opensn
{
namespace lbs
{

class DiffusionDFEMSolver : public LBSSolver
{
protected:
  typedef std::shared_ptr<DiffusionMIPSolver> MIPSolverPtr;

public:
  std::vector<MIPSolverPtr> gs_mip_solvers_;

public:
  explicit DiffusionDFEMSolver(const InputParameters& params);
  ~DiffusionDFEMSolver() override;
  void Initialize() override;
  void InitializeWGSSolvers() override;

public:
  static InputParameters GetInputParameters();
};

} // namespace lbs
} // namespace opensn
