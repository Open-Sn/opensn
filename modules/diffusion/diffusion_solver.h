// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver.h"
#include "modules/diffusion/boundary.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/object_factory.h"

namespace opensn
{

/**
 * Base class for diffusion solvers
 */
class DiffusionSolverBase : public Solver
{
public:
  /** Class for diffusion boundaries */
  class Boundary
  {
  public:
    Boundary();
    Boundary(BoundaryType type, const std::array<double, 3>& values);

    BoundaryType type;
    std::array<double, 3> values;
  };

  explicit DiffusionSolverBase(const std::string& name);
  explicit DiffusionSolverBase(const InputParameters& params);
  ~DiffusionSolverBase() override;

  /**
   * Updates the field functions with the latest data.
   */
  void UpdateFieldFunctions();

protected:
  void InitFieldFunctions();

  using BoundaryInfo = std::pair<BoundaryType, std::vector<double>>;
  using BoundaryPreferences = std::map<std::string, BoundaryInfo>;

  std::shared_ptr<MeshContinuum> grid_ptr_;
  std::shared_ptr<SpatialDiscretization> sdm_ptr_;
  size_t num_local_dofs_;
  size_t num_global_dofs_;
  /// approx solution
  Vec x_;
  /// RHS
  Vec b_;
  /// linear system matrix
  Mat A_;

  std::map<uint64_t, Boundary> boundaries_;
  BoundaryPreferences boundary_preferences_;

public:
  static InputParameters GetInputParameters();
  static InputParameters OptionsBlock();
  static InputParameters BoundaryOptionsBlock();
};

} // namespace opensn
