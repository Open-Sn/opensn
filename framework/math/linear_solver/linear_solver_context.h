// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <petscksp.h>

namespace opensn
{

enum class ResidualScaleType
{
  NONE = 0,
  RHS_NORM = 1,
  RHS_PRECONDITIONED_NORM = 2,
  CUSTOM_SCALE = 3
};

struct LinearSolverContext
{
  virtual int MatrixAction(Mat& matrix, Vec& vector, Vec& action) { return 0; }

  virtual ~LinearSolverContext() = default;
};

struct LinearSystemContext : public LinearSolverContext
{
  double rhs_norm = 0.0;
  double rhs_preconditioned_norm = 0.0;
  double custom_residual_scale = 1.0;
  ResidualScaleType residual_scale_type = ResidualScaleType::NONE;
};

struct LinearEigenContext : public LinearSolverContext
{
  double eigenvalue = 0.0;
};

} // namespace opensn
