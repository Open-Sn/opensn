// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"
#include "framework/parameters/parameter_block.h"
#include <iostream>
#include <utility>

namespace opensn
{
class FieldFunctionGridBased;
class TimeStepper;

class Problem : public Object
{
public:
  explicit Problem(std::string name);
  explicit Problem(const InputParameters& params);
  virtual ~Problem() = default;

  std::string GetName() const;

  std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions();

  const std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions() const;

  /// Initialize function.
  virtual void Initialize();

  /// Generalized query for information supporting varying returns.
  virtual ParameterBlock GetInfo(const ParameterBlock& params) const;

  virtual void SetProperties(const ParameterBlock& params);

  /// PreCheck call to GetInfo.
  ParameterBlock GetInfoWithPreCheck(const ParameterBlock& params) const;

protected:
  std::vector<std::shared_ptr<FieldFunctionGridBased>> field_functions_;

private:
  const std::string name_;

public:
  /// Returns the input parameters.
  static InputParameters GetInputParameters();
};

} // namespace opensn
