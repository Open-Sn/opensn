#pragma once

#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

namespace opensn
{
namespace lbs
{

class DiffusionDFEMSolver : public LBSSolver
{
public:
  std::vector<std::shared_ptr<DiffusionMIPSolver>> gs_mip_solvers_;

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
