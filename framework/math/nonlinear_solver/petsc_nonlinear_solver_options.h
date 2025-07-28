// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"

namespace opensn
{

class PETScNonLinearSolverOptions
{
public:
  explicit PETScNonLinearSolverOptions(const InputParameters& params);

  PETScNonLinearSolverOptions() = default;

  std::string nl_method;
  std::string l_method;
  ParameterBlock pc_options;
  std::string petsc_snes_type;

  double nl_rel_tol;
  double nl_abs_tol;
  double nl_sol_tol;
  int nl_max_its;
  int nl_max_r_evaluations;

  int l_max_failed_iterations;
  double l_rel_tol;
  double l_abs_tol;
  double l_div_tol;
  int l_max_its;
  int l_gmres_restart_intvl;
  double l_gmres_breakdown_tol;

public:
  static InputParameters GetInputParameters();
};

} // namespace opensn
