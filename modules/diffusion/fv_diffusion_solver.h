// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/utils/timer.h"
#include "framework/mesh/mesh.h"
#include "modules/diffusion/diffusion_solver.h"
#include <map>

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;
class ScalarSpatialMaterialFunction;

/// FV diffusion solver
class FVDiffusionSolver : public DiffusionSolverBase
{
public:
  explicit FVDiffusionSolver(const std::string& name, std::shared_ptr<MeshContinuum> grid_ptr);
  explicit FVDiffusionSolver(const InputParameters& params);
  ~FVDiffusionSolver() override;

  void SetBoundaryOptions(const InputParameters& params) override;

  void Initialize() override;
  void Execute() override;

  static InputParameters GetInputParameters();
};

} // namespace opensn
