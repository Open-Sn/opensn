// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include "framework/math/math.h"
#include <utility>

namespace opensn
{

/**
 * An adjoint discrete ordinates solver.
 *
 * This class provides functionality to perform all necessary modifications
 * for an adjoint discrete ordinates solve. This includes transposing the
 * group-to-group transfers, reversing the angles, and evaluating inner
 * products with forward material and point sources for response evaluations.
 *
 * @note In general, distributed sources should be used for volumetric QoIs.
 *       The user is responsible for ensuring that only the appropriate
 *       sources are active in the problem.
 */
class DiscreteOrdinatesAdjointSolver : public DiscreteOrdinatesSolver
{
public:
  explicit DiscreteOrdinatesAdjointSolver(const InputParameters& params);
  explicit DiscreteOrdinatesAdjointSolver(const std::string& solver_name);

  DiscreteOrdinatesAdjointSolver(const DiscreteOrdinatesAdjointSolver&) = delete;
  DiscreteOrdinatesAdjointSolver& operator=(const DiscreteOrdinatesAdjointSolver&) = delete;

  void Initialize() override;
  void Execute() override;

  /// Computes the inner product of the flux and the material source.
  double ComputeInnerProduct();

  /// Exports an importance map in binary format.
  void ExportImportanceMap(const std::string& file_name);

public:
  std::vector<std::vector<double>> flux_moment_buffers_;

public:
  /// Returns the input parameters.
  static InputParameters GetInputParameters();
};

} // namespace opensn
