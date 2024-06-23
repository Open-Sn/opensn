// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/diffusion/diffusion_solver.h"

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;
class ScalarSpatialMaterialFunction;

/// FV diffusion solver
class FVDiffusionSolver : public DiffusionSolverBase
{
public:
  explicit FVDiffusionSolver(const std::string& name);
  explicit FVDiffusionSolver(const InputParameters& params);
  ~FVDiffusionSolver() override;

  void SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);

  void SetOptions(const InputParameters& params);
  void SetBoundaryOptions(const InputParameters& params);

  void Initialize() override;
  void Execute() override;

private:
  std::shared_ptr<ScalarSpatialMaterialFunction> d_coef_function_;
  std::shared_ptr<ScalarSpatialMaterialFunction> sigma_a_function_;
  std::shared_ptr<ScalarSpatialMaterialFunction> q_ext_function_;

public:
  static InputParameters GetInputParameters();
  static InputParameters OptionsBlock();
  static InputParameters BoundaryOptionsBlock();
};

} // namespace opensn
