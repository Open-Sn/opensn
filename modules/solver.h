// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include <iostream>
#include <utility>

namespace opensn
{
class FieldFunctionGridBased;
class TimeStepper;

/**
 * \defgroup SolverBase Base class for all solvers
 * \ingroup doc_PhysicsSolver
 */
class Solver
{
public:
  explicit Solver(std::string name);
  explicit Solver(const InputParameters& params);
  virtual ~Solver() = default;

  std::string GetName() const;

  TimeStepper& GetTimeStepper();
  const TimeStepper& GetTimeStepper() const;

  /// Initialize function.
  virtual void Initialize();

  /// Execution function.
  virtual void Execute();

  /// Step function*/
  virtual void Step();

  /// Advance time values function.
  virtual void Advance();

  /// Generalized query for information supporting varying returns.
  virtual ParameterBlock GetInfo(const ParameterBlock& params) const;

  /// PreCheck call to GetInfo.
  ParameterBlock GetInfoWithPreCheck(const ParameterBlock& params) const;

private:
  const std::string name_;

public:
  /// Returns the input parameters.
  static InputParameters GetInputParameters();
};

} // namespace opensn
