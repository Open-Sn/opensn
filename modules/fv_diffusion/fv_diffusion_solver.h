// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver_base/solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/utils/timer.h"
#include "framework/mesh/mesh.h"
#include "modules/common/diffusion_bndry.h"
#include <map>

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;
class ScalarSpatialMaterialFunction;

namespace diffusion
{

/**
 * FV diffusion solver
 */
class FVSolver : public opensn::Solver
{
public:
  explicit FVSolver(const std::string& name);
  explicit FVSolver(const InputParameters& params);
  ~FVSolver() override;

  void SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);

  void SetOptions(const InputParameters& params);
  void SetBoundaryOptions(const InputParameters& params);

  // void Initialize() override;
  void Initialize() override;
  void Execute() override;

  /**Updates the field functions with the latest data.*/
  void UpdateFieldFunctions();

private:
  typedef std::pair<opensn::diffusion::BoundaryType, std::vector<double>> BoundaryInfo;
  typedef std::map<std::string, BoundaryInfo> BoundaryPreferences;

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

  std::shared_ptr<ScalarSpatialMaterialFunction> d_coef_function_;
  std::shared_ptr<ScalarSpatialMaterialFunction> sigma_a_function_;
  std::shared_ptr<ScalarSpatialMaterialFunction> q_ext_function_;

public:
  static InputParameters GetInputParameters();
  static InputParameters OptionsBlock();
  static InputParameters BoundaryOptionsBlock();
};

} // namespace diffusion
} // namespace opensn
