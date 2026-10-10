// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include <memory>
#include <string>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

/**
 * Base class for solvers that drive a DiscreteOrdinatesProblem.
 *
 * Solution and time state remain in the shared problem. A solver restores its temporary settings
 * before returning. Validation is collective.
 */
class DiscreteOrdinatesSolver : public Solver
{
protected:
  explicit DiscreteOrdinatesSolver(const InputParameters& params);

  /// Requires a built problem, validates its configuration, then checks solver requirements.
  void ValidateState() const final;

  /**
   * Appends one message for each solver-specific requirement the problem violates.
   *
   * Problem-wide rules belong in DiscreteOrdinatesProblem::CheckConfigurationErrors. Requirements
   * based on rank-local state reduce it before appending an error.
   */
  virtual void CheckRequirements(std::vector<std::string>& errors) const = 0;

  /// Appends the error for a solver that requires steady-state mode.
  void RequireSteadyState(std::vector<std::string>& errors) const;

  /// Appends an error for each external source (volumetric, point, or incoming-flux boundary).
  void RequireNoExternalSources(std::vector<std::string>& errors) const;

  /// Appends an error if no cell has fissionable cross sections. Collective.
  void RequireFissionableMaterial(std::vector<std::string>& errors) const;

  /// Appends the requirements common to k-eigenvalue solvers. Collective.
  void CheckKEigenRequirements(std::vector<std::string>& errors) const;

  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;

public:
  static InputParameters GetInputParameters();
};

} // namespace opensn
