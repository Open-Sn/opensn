// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"

namespace opensn
{

class NonLinearSolverOptions : public Object
{
public:
  static InputParameters GetInputParameters();
  explicit NonLinearSolverOptions(const InputParameters& params);
  NonLinearSolverOptions() = default;

  std::string nl_method = "JFNK";
  std::string l_method = "gmres";

  ParameterBlock pc_options;

  std::string petsc_snes_type = "newtonls";

  double nl_rel_tol = 1.0e-8;
  double nl_abs_tol = 1.0e-8;
  double nl_sol_tol = 1.0e-50;
  int nl_max_its = 50;
  int nl_max_r_evaluations = -1;
  int l_max_failed_iterations = 1000;
  double l_rel_tol = 1.0e-8;
  double l_abs_tol = 1.0e-8;
  double l_div_tol = 1.0e6;
  int l_max_its = 100;
  int l_gmres_restart_intvl = 30;
  double l_gmres_breakdown_tol = 1.0e6;
};

} // namespace opensn
