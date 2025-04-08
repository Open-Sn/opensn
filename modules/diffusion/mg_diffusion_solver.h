// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/diffusion/boundary.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/physics/solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/utils/timer.h"
#include "framework/math/vector.h"
#include <map>
#include <set>

namespace opensn
{
class MeshContinuum;
class SpatialDiscretization;

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

/// Multi-group diffusion solver
class MGDiffusionSolver : public Solver
{
public:
  /** Multigroup diffusion boundary */
  class MGBoundary
  {
  public:
    MGBoundary();
    MGBoundary(BoundaryType type, const std::array<std::vector<double>, 3>& mg_values);

    BoundaryType type;
    std::array<std::vector<double>, 3> mg_values;
  };

  explicit MGDiffusionSolver(const std::string& name, std::shared_ptr<MeshContinuum> grid_ptr);
  explicit MGDiffusionSolver(const InputParameters& params);
  ~MGDiffusionSolver() override;

  void SetOptions(const InputParameters& params);
  void SetBoundaryOptions(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

  /// Updates the field functions with the latest data.
  void UpdateFieldFunctions();

  const std::vector<std::shared_ptr<FieldFunctionGridBased>>& GetFieldFunctions() const;

private:
  void InitializeMaterials(std::set<int>& material_ids);
  void SetBCs(const std::vector<uint64_t>& global_unique_bndry_ids);
  void AssembleAbext();
  void ComputeTwoGridParams();
  void ComputeTwoGridVolumeFractions();
  void AssembleRhs(unsigned int g, int64_t iverbose);
  void AssembleRhsTwoGrid(int64_t iverbose);
  void SolveOneGroupProblem(unsigned int g, int64_t iverbose);
  void UpdateFluxWithTwoGrid(int64_t iverbose);

  using BoundaryInfo = std::pair<BoundaryType, std::array<std::vector<double>, 3>>;
  using BoundaryPreferences = std::map<std::string, BoundaryInfo>;

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
  std::map<int, std::shared_ptr<MultiGroupXS>> block_id_to_xs_map_;

  std::map<int, TwoGridCollapsedInfo> map_mat_id_2_tginfo_;
  //  std::map<int, Multigroup_D_and_sigR> map_mat_id_2_tgXS;

  std::vector<std::shared_ptr<FieldFunctionGridBased>> field_functions_;

public:
  static InputParameters GetInputParameters();
  static InputParameters GetOptionsBlock();
  static InputParameters GetBoundaryOptionsBlock();
  static InputParameters GetXSMapEntryBlock();
};

} // namespace opensn
