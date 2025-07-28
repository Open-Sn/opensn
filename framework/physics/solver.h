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
  /// Returns the input parameters.
  static InputParameters GetInputParameters();
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

  /**
   * \addtogroup SolverBase
   *
   * \section Properties Properties that can be set via
   * `SetProperties`
   * \copydoc opensn::Solver::SetProperties
   *
   * Base solver settable properties:
   * - `dt`, Timestep size
   * - `time`, Current time
   */
  virtual void SetProperties(const ParameterBlock& params);

  /// PreCheck call to GetInfo.
  ParameterBlock GetInfoWithPreCheck(const ParameterBlock& params) const;

protected:
  std::shared_ptr<TimeStepper> timestepper_ = nullptr;

private:
  const std::string name_;

private:
  static std::shared_ptr<TimeStepper> InitTimeStepper(const InputParameters& params);
};

} // namespace opensn
