// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include "framework/object.h"

namespace opensn
{
class LBSProblem;
class LBSGroupset;
class PowerIterationKEigenSolver;

/**
 * Base class for LBS acceleration methods.
 */
class LBSAcceleration : public Object
{
public:
  static InputParameters GetInputParameters();

  explicit LBSAcceleration(const InputParameters& params);

  virtual ~LBSAcceleration() = default;

  /**
   * Public initialize method, to be called by the power iteration solver
   * using this acceleration scheme.
   *
   * Calls the internal Initialize() method.
   */
  void Initialize(PowerIterationKEigenSolver& solver);

  /**
   * Initialization method for derived classes. Must be overridden.
   *
   * Called by the other Initialize() method that takes a solver.
   */
  virtual void Initialize() = 0;

  /**
   * Pre-execute method for derived classes. Must be overridden.
   *
   * Called by the owning solver at the beginning of its Execute().
   */
  virtual void PreExecute() = 0;

  /**
   * Pre-power iteration method for derived classes. Must be overridden.
   *
   * Called by the owning solver at the beginning of the power iteration loop
   * after the fission source is set and the local moments are scaled.
   */
  virtual void PrePowerIteration() = 0;

  /**
   * Post-power iteration method for derived classes. Must be overridden.
   *
   * Called by the owning solver at the end of the power iteration loop
   * after transport is solved.
   *
   * Returns the k-eigenvalue from the acceleration solve.
   */
  virtual double PostPowerIteration() = 0;

  const std::string& GetName() const { return name_; }

protected:
  /// The associated LBS problem
  LBSProblem& lbs_problem_;

  /// Absolute residual tolerance from parameters
  const double l_abs_tol_;
  /// Maximum allowable iterations from parameters
  const int max_iters_;
  /// Verbosity flag from parameters
  const bool verbose_;
  /// PETSc options from parameters
  const std::string petsc_options_;

  /// Groupsets from the LBSProblem
  std::vector<LBSGroupset>& groupsets_;
  /// Source moments vector from the LBSProblem
  std::vector<double>& q_moments_local_;
  /// Last updated flux vector from the LBSProblem
  std::vector<double>& phi_old_local_;
  /// Newest updated flux vector from the LBSProblem
  std::vector<double>& phi_new_local_;

  /// Associated PI solver, filled during Initialize()
  PowerIterationKEigenSolver* solver_;

private:
  const std::string name_;
};
} // namespace opensn
