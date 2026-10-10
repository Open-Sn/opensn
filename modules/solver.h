// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include <iostream>
#include <utility>

namespace opensn
{
class FieldFunctionGridBased;

/**
 * Base class for solver lifecycle management.
 *
 * Public lifecycle functions validate state before calling protected derived-class hooks.
 * Execute() and Advance() also require successful initialization.
 */
class Solver
{
public:
  explicit Solver(std::string name);
  explicit Solver(const InputParameters& params);
  virtual ~Solver() = default;

  std::string GetName() const;

  /**
   * Validates and initializes the solver.
   *
   * The solver is marked initialized only after InitializeSolver() and post-initialization
   * validation succeed. Calling this again clears the previous initialized state first.
   */
  void Initialize();

  /// Validates the current state and executes the solver. Requires successful Initialize().
  void Execute();

  /// Validates the current state and advances one step. Requires successful Initialize().
  void Advance();

  /// Returns true once Initialize() has completed.
  bool IsInitialized() const { return initialized_; }

  /// Generalized query for information supporting varying returns.
  virtual ParameterBlock GetInfo(const ParameterBlock& params) const;

  /// PreCheck call to GetInfo.
  ParameterBlock GetInfoWithPreCheck(const ParameterBlock& params) const;

  bool IsBalanceEnabled() const { return compute_balance_; }

protected:
  /**
   * Validates the solver's options and the object it drives.
   *
   * Called before each lifecycle hook and again after InitializeSolver(). Implementations throw
   * std::invalid_argument for unsupported configurations. MPI implementations are collective.
   */
  virtual void ValidateState() const = 0;

  /// Performs derived-class initialization after ValidateState() succeeds.
  virtual void InitializeSolver() = 0;
  /// Performs derived-class execution after initialization and validation succeed.
  virtual void ExecuteSolver() = 0;
  /// Performs one derived-class step. The default logs that Advance() is unsupported.
  virtual void AdvanceSolver();

private:
  const std::string name_;
  bool compute_balance_ = false;
  bool initialized_ = false;

public:
  /// Returns the input parameters.
  static InputParameters GetInputParameters();
};

} // namespace opensn
