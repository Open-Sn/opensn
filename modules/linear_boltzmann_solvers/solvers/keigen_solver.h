// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver.h"

namespace opensn
{

class DiscreteOrdinatesProblem;
class FieldFunctionGridBased;

class KEigenSolver : public Solver
{
public:
  explicit KEigenSolver(const InputParameters& params);

  /// Returns the power generation field function, if enabled.
  std::shared_ptr<FieldFunctionGridBased> GetPowerFieldFunction() const;

  void UpdatePowerFieldFunctions();

protected:
  void InitializePowerFieldFunction();

  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;

  std::shared_ptr<FieldFunctionGridBased> power_field_func_;

public:
  static InputParameters GetInputParameters();
};

} // namespace opensn
