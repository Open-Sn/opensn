// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include <iostream>
#include <utility>

namespace opensn
{
class FieldFunctionGridBased;

class Solver
{
public:
  explicit Solver(std::string name);
  explicit Solver(const InputParameters& params);
  virtual ~Solver() = default;

  std::string GetName() const;

  /// Initialize function.
  virtual void Initialize();

  /// Execution function.
  virtual void Execute();

  /// Advance time values function.
  virtual void Advance();

  /// Generalized query for information supporting varying returns.
  virtual ParameterBlock GetInfo(const ParameterBlock& params) const;

  /// PreCheck call to GetInfo.
  ParameterBlock GetInfoWithPreCheck(const ParameterBlock& params) const;

  bool IsBalanceEnabled() const { return compute_balance_; }

private:
  const std::string name_;
  bool compute_balance_ = false;

public:
  /// Returns the input parameters.
  static InputParameters GetInputParameters();
};

} // namespace opensn
