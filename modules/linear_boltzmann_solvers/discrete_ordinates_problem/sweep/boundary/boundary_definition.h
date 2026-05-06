// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/function.h"
#include "framework/parameters/input_parameters.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include <limits>

namespace opensn
{

/**
 * Parsed boundary-condition specification for sweep boundary construction.
 *
 * This structure stores the user-supplied boundary type and any type-specific data needed later by
 * the sweep boundary factory. It intentionally describes the boundary condition without owning a
 * concrete SweepBoundary object.
 */
struct BoundaryDefinition
{
  /// Construct a vacuum boundary definition.
  BoundaryDefinition() = default;
  /**
   * Construct a boundary definition from input parameters.
   *
   * The constructor validates the requested boundary type and required/forbidden type-specific
   * parameters. Isotropic boundaries store group strengths, arbitrary boundaries store an
   * AngularFluxFunction, and vacuum/reflecting boundaries do not store additional payload.
   * \param params Boundary option parameters containing at least `name` and `type`.
   * \param num_groups Number of solver groups used to validate isotropic strengths.
   */
  BoundaryDefinition(const InputParameters& params, unsigned int num_groups);

  /// Boundary condition type.
  LBSBoundaryType type = LBSBoundaryType::VACUUM;
  /// Group-wise incoming angular flux for isotropic boundaries.
  std::vector<double> group_strength;
  // Start time.
  double start_time = -std::numeric_limits<double>::infinity();
  // End time.
  double end_time = std::numeric_limits<double>::infinity();
  /// User-supplied angular flux callback for time-dependent arbitrary boundaries.
  std::shared_ptr<AngularFluxTimeFunction> time_angular_flux_function;
};

} // namespace opensn
