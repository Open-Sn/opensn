// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/problems/diffusion/diffusion_solver.h"
#include "framework/mesh/mesh.h"
#include "physics/solvers/solver.h"
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

  void SetBoundaryOptions(const InputParameters& params) override;

  void Initialize() override;
  void Execute() override;

  static InputParameters GetInputParameters();
  static std::shared_ptr<CFEMDiffusionSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
