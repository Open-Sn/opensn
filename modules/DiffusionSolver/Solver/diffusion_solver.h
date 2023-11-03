#pragma once

#include "mesh/Cell/cell.h"

#include "DiffusionSolver/Boundaries/chi_diffusion_bndry.h"
#include "DiffusionSolver/UnitIntegralContainer.h"

#include "DiffusionSolver/chi_diffusion.h"
#include "physics/SolverBase/chi_solver.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"
#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"

#include "mesh/VolumeMesher/chi_volumemesher.h"

#include "utils/chi_timer.h"

#include <petscksp.h>

#define DIFFUSION_MATERIALS_REGULAR 10
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTR 11
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF 12
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JPART 13
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JFULL 14

//######################################################### Class def
/**Solver for the general diffusion problem.

 <img src="DiffusionMatProp.png" style="width:500px">
 * */
class chi_diffusion::Solver : public chi_physics::Solver
{
private:
  chi::Timer t_assembly_;
  chi::Timer t_solve_;

  double time_assembly_ = 0.0;
  double time_solve_ = 0.0;
  bool verbose_info_ = true;

public:
  typedef unsigned int uint;
  typedef std::pair<BoundaryType, std::vector<double>> BoundaryInfo;
  typedef std::map<std::string, BoundaryInfo> BoundaryPreferences;

public:
  BoundaryPreferences boundary_preferences_;
  std::map<uint64_t, chi_diffusion::Boundary*> boundaries_;
  chi_mesh::MeshContinuumPtr grid_ptr_ = nullptr;

  std::shared_ptr<chi_math::SpatialDiscretization> discretization_;

  chi_math::UnknownManager unknown_manager_;

  int material_mode_ = DIFFUSION_MATERIALS_REGULAR;

  bool common_items_initialized_ = false;

  Vec x_ = nullptr;   // approx solution
  Vec b_ = nullptr;   // RHS
  Mat A_ = nullptr;   // linear system matrix
  KSP ksp_ = nullptr; // linear solver context
  PC pc_ = nullptr;   // preconditioner context

  PetscReal norm_ = 0.0;    /* norm of solution error */
  PetscErrorCode ierr_ = 0; // General error code

  size_t local_dof_count_ = 0;
  size_t global_dof_count_ = 0;

  std::vector<double> pwld_phi_local_;

  int gi_ = 0;
  int G_ = 1;
  std::string options_string_;

  std::map<uint64_t, UnitIntegralContainer> unit_integrals_;

public:
  // 00
  Solver(const Solver&) = delete;
  Solver& operator=(const Solver&) = delete;

  explicit Solver(const std::string& in_solver_name);
  virtual ~Solver();

  // 01 General
  void GetMaterialProperties(const chi_mesh::Cell& cell,
                             int cell_dofs,
                             std::vector<double>& diffCoeff,
                             std::vector<double>& sourceQ,
                             std::vector<double>& sigmaa,
                             int group = 0,
                             int moment = 0);
  // 01a
  void InitializeCommonItems();

  // 01b
  void Initialize() override { Initialize(true); }
  int Initialize(bool verbose);

  // 02a
  void Execute() override { ExecuteS(); }
  int ExecuteS(bool suppress_assembly = false, bool suppress_solve = false);

  void CFEM_Assemble_A_and_b(chi_mesh::Cell& cell, int group = 0);

  // 02c_c
  void PWLD_Assemble_A_and_b(const chi_mesh::Cell& cell, int component = 0);
  void PWLD_Assemble_b(const chi_mesh::Cell& cell, int component = 0);

  // 02e_c
  //  void PWLD_Assemble_A_and_b_GAGG(const chi_mesh::Cell& cell);
  //  void PWLD_Assemble_b_GAGG(const chi_mesh::Cell& cell);

  // 03b
  double HPerpendicular(const chi_mesh::Cell& cell,
                        const UnitIntegralContainer& fe_intgrl_values,
                        unsigned int f);

  static uint64_t MapCellLocalNodeIDFromGlobalID(const chi_mesh::Cell& cell,
                                                 uint64_t node_global_id);

  static unsigned int
  MapCellFace(const chi_mesh::Cell& cur_cell, const chi_mesh::Cell& adj_cell, unsigned int f);

  void UpdateFieldFunctions();
};
