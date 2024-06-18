// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/diffusion/boundary.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/isotropic_multigroup_source.h"
#include "framework/physics/solver_base/solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/utils/timer.h"
#include <map>
#include <set>

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;

namespace diffusion
{

struct KSPAppContext
{
  PetscBool verbose = PETSC_FALSE;
};

struct TwoGridCollapsedInfo
{
  double collapsed_D = 0.0;
  double collapsed_sig_a = 0.0;
  std::vector<double> spectrum;
};

// struct Multigroup_D_and_sigR
//{
//   std::vector<double> Dg;
//   std::vector<double> sigR;
// };

/**
 * Multi-group diffusion solver
 */
class MGSolver : public opensn::Solver
{
public:
  explicit MGSolver(const std::string& name);
  explicit MGSolver(const InputParameters& params);
  ~MGSolver() override;

  void SetOptions(const InputParameters& params);
  void SetBoundaryOptions(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

  /**
   * Updates the field functions with the latest data.
   */
  void UpdateFieldFunctions();

private:
  void InitializeMaterials(std::set<int>& material_ids);
  void SetBCs(const std::vector<uint64_t>& globl_unique_bndry_ids);
  void AssembleAbext();
  void ComputeTwoGridParams();
  void ComputeTwoGridVolumeFractions();
  void AssembleRhs(unsigned int g, int64_t iverbose);
  void AssembleRhsTwoGrid(int64_t iverbose);
  void SolveOneGroupProblem(unsigned int g, int64_t iverbose);
  void UpdateFluxWithTwoGrid(int64_t iverbose);

  typedef std::pair<opensn::diffusion::BoundaryType, std::array<std::vector<double>, 3>>
    BoundaryInfo;

  typedef std::map<std::string, BoundaryInfo> BoundaryPreferences;

  std::shared_ptr<MeshContinuum> grid_ptr_;

  std::shared_ptr<SpatialDiscretization> sdm_ptr_;

  uint num_groups_;
  uint last_fast_group_;
  bool do_two_grid_;

  size_t num_local_dofs_;
  size_t num_global_dofs_;

  /// linear system matrix for each group
  std::vector<Mat> A_;
  /// external source vector for each group
  std::vector<Vec> bext_;
  /// solution vector for each group
  std::vector<Vec> x_;
  /// vector of old fluxes
  std::vector<Vec> x_old_;

  /// error vector for thermal fluxes
  Vec thermal_dphi_;
  /// actual rhs vector for the linear system A[g] x[g] = b
  Vec b_;

  PETScSolverSetup petsc_solver_;
  KSPAppContext my_app_context_;

  std::vector<std::vector<double>> VF_;

  BoundaryPreferences boundary_preferences_;
  std::map<uint64_t, MGBoundary> boundaries_;

protected:
  std::map<int, std::shared_ptr<MultiGroupXS>> matid_to_xs_map;

  std::map<int, std::shared_ptr<IsotropicMultiGroupSource>> matid_to_src_map;

  std::map<int, TwoGridCollapsedInfo> map_mat_id_2_tginfo;
  //  std::map<int, Multigroup_D_and_sigR> map_mat_id_2_tgXS;

public:
  static InputParameters GetInputParameters();
  static InputParameters OptionsBlock();
  static InputParameters BoundaryOptionsBlock();
};

} // namespace diffusion
} // namespace opensn
