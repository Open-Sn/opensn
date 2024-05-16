// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver_base/solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"

#include "modules/cfem_diffusion/cfem_diffusion_bndry.h"
#include "framework/utils/timer.h"
#include <map>

#include "framework/mesh/mesh.h"

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;
class ScalarSpatialMaterialFunction;

namespace diffusion
{

/**
 * CFEM diffusion solver
 *
 */
class CFEMSolver : public opensn::Solver
{
public:
  explicit CFEMSolver(const std::string& name);
  explicit CFEMSolver(const InputParameters& params);
  ~CFEMSolver() override;

  void SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);

  void SetOptions(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

  /**
   * Updates the field functions with the latest data.
   */
  void UpdateFieldFunctions();

public:
  typedef std::pair<opensn::diffusion::BoundaryType, std::vector<double>> BoundaryInfo;
  typedef std::map<std::string, BoundaryInfo> BoundaryPreferences;
  BoundaryPreferences boundary_preferences;

private:
  std::shared_ptr<MeshContinuum> grid_ptr_ = nullptr;

  std::shared_ptr<SpatialDiscretization> sdm_ptr_ = nullptr;

  size_t num_local_dofs_ = 0;
  size_t num_globl_dofs_ = 0;

  /// approx solution
  Vec x_ = nullptr;
  /// RHS
  Vec b_ = nullptr;
  /// linear system matrix
  Mat A_ = nullptr;

  std::map<uint64_t, Boundary> boundaries_;

  std::shared_ptr<ScalarSpatialMaterialFunction> d_coef_function_;
  std::shared_ptr<ScalarSpatialMaterialFunction> sigma_a_function_;
  std::shared_ptr<ScalarSpatialMaterialFunction> q_ext_function_;

public:
  static InputParameters GetInputParameters();
  static InputParameters OptionsBlock();
};

} // namespace diffusion
} // namespace opensn
