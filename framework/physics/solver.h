// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"
#include "framework/physics/basic_options.h"
#include "framework/parameters/parameter_block.h"
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
class Solver : public Object
{
public:
  explicit Solver(std::string name);
  Solver(std::string name, std::initializer_list<BasicOption> options);
  explicit Solver(const InputParameters& params);

  virtual ~Solver() = default;

  std::string Name() const;

  BasicOptions& GetBasicOptions();
  const BasicOptions& GetBasicOptions() const;

  std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions();
  const std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions() const;

  TimeStepper& GetTimeStepper();
  const TimeStepper& GetTimeStepper() const;

  virtual void Initialize();
  virtual void Execute();
  /// Solve the current time step.
  virtual void Step();
  /// Move the solver to the start of next time step.
  virtual void Advance();

  /// Generalized query for information supporting varying returns.
  virtual ParameterBlock GetInfo(const ParameterBlock& params) const;

  /**
   * \addtogroup SolverBase
   *
   * \section Properties Properties that can be set
   * The following properties can be set via the lua call
   * `SolverSetProperties`
   * \copydoc opensn::Solver::SetProperties
   *
   * Base solver settable properties:
   * - `dt`, Time step size
   * - `time`, Current time
   */
  virtual void SetProperties(const ParameterBlock& params);

  /// Pre-check call to GetInfo.
  ParameterBlock GetInfoWithPreCheck(const ParameterBlock& params) const;

protected:
  void SetAuxiliaryFieldFunction();

protected:
  BasicOptions basic_options_;
  std::vector<std::shared_ptr<FieldFunctionGridBased>> field_functions_;
  std::map<std::string, std::shared_ptr<FieldFunctionGridBased>> aux_field_functions_;
  std::shared_ptr<TimeStepper> timestepper_ = nullptr;

private:
  const std::string name_;
  const std::vector<std::string> auxvars_;

public:
  static InputParameters GetInputParameters();

private:
  static std::shared_ptr<TimeStepper> InitializeTimeStepper(const InputParameters& params);
};

} // namespace opensn
