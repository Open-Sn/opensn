// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/diffusion/diffusion_solver.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/mesh/mesh.h"
#include <vector>

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;
class ScalarSpatialMaterialFunction;

/// DFEM diffusion solver
class DFEMDiffusionSolver : public DiffusionSolverBase
{
public:
  explicit DFEMDiffusionSolver(const std::string& name);
  explicit DFEMDiffusionSolver(const InputParameters& params);
  ~DFEMDiffusionSolver() override;

  void SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);

  void SetOptions(const InputParameters& params);
  void SetBoundaryOptions(const InputParameters& params);

  void Initialize() override;
  void Execute() override;

private:
  /**
   * Still searching for a reference for this.
   *
   * For Polygons:
   * Defined from paper
   * Turcksin B, Ragusa J, "Discontinuous diffusion synthetic acceleration for S_n transport on 2D
   * arbitrary polygonal meshes", Journal of Computational Physics 274, pg 356-369, 2014.
   * Nv =  Number of vertices. If Nv <= 4 then the perimeter parameter should be replaced by edge
   * length.
   */
  double HPerpendicular(const Cell& cell, unsigned int f);

  /// Maps a face, in a discontinuous sense, using the spatial discretization.
  int MapFaceNodeDisc(const Cell& cur_cell,
                      const Cell& adj_cell,
                      const std::vector<Vector3>& cc_node_locs,
                      const std::vector<Vector3>& ac_node_locs,
                      size_t ccf,
                      size_t acf,
                      size_t ccfi,
                      double epsilon = 1.0e-12);

  std::vector<double> field_;

  std::shared_ptr<ScalarSpatialMaterialFunction> d_coef_function_;
  std::shared_ptr<ScalarSpatialMaterialFunction> sigma_a_function_;
  std::shared_ptr<ScalarSpatialMaterialFunction> q_ext_function_;

public:
  static InputParameters GetInputParameters();
  static InputParameters OptionsBlock();
  static InputParameters BoundaryOptionsBlock();
  static std::shared_ptr<DFEMDiffusionSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
