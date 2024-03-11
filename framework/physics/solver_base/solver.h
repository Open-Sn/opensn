#pragma once

#include "framework/object.h"
#include "framework/physics/physics_namespace.h"

#include "framework/physics/basic_options/basic_options.h"
#include "framework/parameters/parameter_block.h"

#include <iostream>
#include <utility>

namespace opensn
{
class FieldFunctionGridBased;
class TimeStepper;

/**\defgroup SolverBase Base class for all solvers
 * \ingroup doc_PhysicsSolver*/
class Solver : public Object
{
public:
  /**Returns the input parameters.*/
  static InputParameters GetInputParameters();
  explicit Solver(std::string name);
  Solver(std::string name, std::initializer_list<BasicOption> options);
  explicit Solver(const InputParameters& params);
  virtual ~Solver() = default;

  std::string TextName() const;

  BasicOptions& GetBasicOptions();
  const BasicOptions& GetBasicOptions() const;

  std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions();

  const std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions() const;

  TimeStepper& GetTimeStepper();
  const TimeStepper& GetTimeStepper() const;

  /**Initialize function.*/
  virtual void Initialize();
  /**Execution function.*/
  virtual void Execute();
  /**Step function*/
  virtual void Step();
  /**Advance time values function.*/
  virtual void Advance();

  /**Generalized query for information supporting varying returns.*/
  virtual ParameterBlock GetInfo(const ParameterBlock& params) const;
  /**\addtogroup SolverBase
   *
   * \section Properties Properties that can be set
   * The following properties can be set via the lua call
   * `SolverSetProperties`
   * \copydoc opensn::Solver::SetProperties
   *
   * Base solver settable properties:
   * - `dt`, Timestep size
   * - `time`, Current time
   */
  virtual void SetProperties(const ParameterBlock& params);
  /**PreCheck call to GetInfo.*/
  ParameterBlock GetInfoWithPreCheck(const ParameterBlock& params) const;

protected:
  BasicOptions basic_options_;
  std::vector<std::shared_ptr<FieldFunctionGridBased>> field_functions_;
  std::shared_ptr<TimeStepper> timestepper_ = nullptr;

private:
  static std::shared_ptr<TimeStepper> InitTimeStepper(const InputParameters& params);
  const std::string text_name_;
};

} // namespace opensn
