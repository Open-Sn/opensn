#pragma once

#include "opensn/framework/ChiObject.h"
#include "opensn/framework/physics/chi_physics_namespace.h"

#include "opensn/framework/physics/BasicOptions/basic_options.h"
#include "opensn/framework/parameters/parameter_block.h"

#include <iostream>
#include <utility>

namespace chi_physics
{
class FieldFunctionGridBased;
class TimeStepper;

// ######################################################### Solver parent class
/**\defgroup SolverBase Base class for all solvers
 * \ingroup doc_PhysicsSolver*/
class Solver : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();
  explicit Solver(std::string in_text_name);
  Solver(std::string in_text_name, std::initializer_list<BasicOption> in_options);
  explicit Solver(const chi::InputParameters& params);
  virtual ~Solver() = default;

  std::string TextName() const;

  BasicOptions& GetBasicOptions();
  const BasicOptions& GetBasicOptions() const;

  std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions();

  const std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions() const;

  TimeStepper& GetTimeStepper();
  const TimeStepper& GetTimeStepper() const;

  virtual void Initialize();
  virtual void Execute();
  virtual void Step();
  virtual void Advance();

  /**Generalized query for information supporting varying returns.*/
  virtual chi::ParameterBlock GetInfo(const chi::ParameterBlock& params) const;
  virtual void SetProperties(const chi::ParameterBlock& params);
  /**PreCheck call to GetInfo.*/
  chi::ParameterBlock GetInfoWithPreCheck(const chi::ParameterBlock& params) const;

protected:
  BasicOptions basic_options_;
  std::vector<std::shared_ptr<FieldFunctionGridBased>> field_functions_;
  std::shared_ptr<TimeStepper> timestepper_ = nullptr;

private:
  static std::shared_ptr<TimeStepper> InitTimeStepper(const chi::InputParameters& params);
  const std::string text_name_;
};

} // namespace chi_physics
