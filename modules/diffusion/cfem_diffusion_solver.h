// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/diffusion/diffusion_solver.h"
#include "framework/mesh/mesh.h"
#include "framework/physics/solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/utils/timer.h"
#include <map>

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;
class ScalarSpatialMaterialFunction;

/**
 * CFEM diffusion solver
 */
class CFEMDiffusionSolver : public DiffusionSolverBase
{
public:
  explicit CFEMDiffusionSolver(const std::string& name, std::shared_ptr<MeshContinuum> grid_ptr);
  explicit CFEMDiffusionSolver(const InputParameters& params);
  ~CFEMDiffusionSolver() override;

  void SetDCoefFunction(ScalarSpatialMaterialFunction* function);
  void SetQExtFunction(ScalarSpatialMaterialFunction* function);
  void SetSigmaAFunction(ScalarSpatialMaterialFunction* function);

  void SetOptions(const InputParameters& params);
  void SetBoundaryOptions(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

private:
  ScalarSpatialMaterialFunction* d_coef_function_;
  ScalarSpatialMaterialFunction* sigma_a_function_;
  ScalarSpatialMaterialFunction* q_ext_function_;

public:
  static InputParameters GetInputParameters();
  static InputParameters GetOptionsBlock();
  static InputParameters GetBoundaryOptionsBlock();
  static std::shared_ptr<CFEMDiffusionSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
