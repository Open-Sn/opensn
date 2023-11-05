#pragma once

#include "opensn/framework/ChiObject.h"

namespace chi_math
{

class NonLinearSolverOptions : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();
  explicit NonLinearSolverOptions(const chi::InputParameters& params);
  NonLinearSolverOptions() = default;

  std::string nl_method_ = "JFNK";
  std::string l_method_ = "gmres";

  chi::ParameterBlock pc_options_;

  std::string petsc_snes_type_ = "newtonls";

  double nl_rel_tol_ = 1.0e-8;
  double nl_abs_tol_ = 1.0e-8;
  double nl_sol_tol_ = 1.0e-50;
  int nl_max_its_ = 50;
  int nl_max_r_evaluations_ = -1;
  int l_max_failed_iterations_ = 1000;
  double l_rel_tol_ = 1.0e-8;
  double l_abs_tol_ = 1.0e-8;
  double l_div_tol_ = 1.0e6;
  int l_max_its_ = 100;
  int l_gmres_restart_intvl_ = 30;
  double l_gmres_breakdown_tol_ = 1.0e6;
};

} // namespace chi_math
