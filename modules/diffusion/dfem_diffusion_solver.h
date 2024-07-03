// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver_base/solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "modules/diffusion/boundary.h"
#include "framework/utils/timer.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/mesh/mesh.h"
#include <map>

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;
class ScalarSpatialMaterialFunction;

namespace diffusion
{

/**
 * DFEM diffusion solver
 */
class DFEMSolver : public opensn::Solver
{
public:
  explicit DFEMSolver(const std::string& name);
  explicit DFEMSolver(const InputParameters& params);
  ~DFEMSolver() override;

  void SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);
  void SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function);

  void SetOptions(const InputParameters& params);
  void SetBoundaryOptions(const InputParameters& params);

  void Initialize() override;
  void Execute() override;

  /**
   * Updates the field functions with the latest data.
   */
  void UpdateFieldFunctions();

private:
  /**Still searching for a reference for this.
   *
   * For Polygons:
   * Defined from paper  \n
   * Turcksin B, Ragusa J, "Discontinuous diffusion synthetic acceleration
   * for S_n transport on 2D arbitrary polygonal meshes", Journal of
   * Computational Physics 274, pg 356-369, 2014.\n
   * \n
   * Nv = Number of vertices. If Nv <= 4 then the perimeter parameter
   * should be replaced by edge length.*/
  double HPerpendicular(const Cell& cell, unsigned int f);

  /**
   * Maps a face, in a discontinuous sense, using the spatial discretization.
   */
  int MapFaceNodeDisc(const Cell& cur_cell,
                      const Cell& adj_cell,
                      const std::vector<Vector3>& cc_node_locs,
                      const std::vector<Vector3>& ac_node_locs,
                      size_t ccf,
                      size_t acf,
                      size_t ccfi,
                      double epsilon = 1.0e-12);

  typedef std::pair<opensn::diffusion::BoundaryType, std::vector<double>> BoundaryInfo;
  typedef std::map<std::string, BoundaryInfo> BoundaryPreferences;

  std::shared_ptr<MeshContinuum> grid_ptr_ = nullptr;

  std::shared_ptr<SpatialDiscretization> sdm_ptr_ = nullptr;

  size_t num_local_dofs_ = 0;
  size_t num_globl_dofs_ = 0;

  std::vector<double> field_;

  /// approx solution
  Vec x_ = nullptr;
  /// RHS
  Vec b_ = nullptr;
  /// linear system matrix
  Mat A_ = nullptr;

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
