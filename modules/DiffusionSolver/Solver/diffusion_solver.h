#pragma once

#include "opensn/framework/mesh/Cell/cell.h"
#include "opensn/modules/DiffusionSolver/Boundaries/chi_diffusion_bndry.h"
#include "opensn/modules/DiffusionSolver/UnitIntegralContainer.h"
#include "opensn/modules/DiffusionSolver/chi_diffusion.h"
#include "opensn/framework/physics/SolverBase/chi_solver.h"
#include "opensn/framework/math/SpatialDiscretization/SpatialDiscretization.h"
#include "opensn/framework/math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"
#include "opensn/framework/math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"
#include "opensn/framework/mesh/VolumeMesher/chi_volumemesher.h"
#include "opensn/framework/utils/chi_timer.h"
#include <petscksp.h>

#define DIFFUSION_MATERIALS_REGULAR 10
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTR 11
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF 12
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JPART 13
#define DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JFULL 14

namespace chi_diffusion
{

/**
 * Solver for the general diffusion problem.
 *
 * <img src="DiffusionMatProp.png" style="width:500px">
 *
 */
class Solver : public chi_physics::Solver
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
  std::map<uint64_t, Boundary*> boundaries_;
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
  Solver(const Solver&) = delete;
  Solver& operator=(const Solver&) = delete;

  /**
   * \defgroup LuaDiffusionBasicOptions Basic Options
   * \ingroup LuaDiffusion
   *
   * Option name           | Type   | Default Value | Description
   * ----------------------|------- |---------------|------------
   * "discretization_method" | string | "None"        | Spatial discretization
   * method. "PWLC", "PWLD_MIP". "max_iters"             | int    | 500           |
   * Maximum iterations for solver convergence. "residual_tolerance"    | float
   * | 1.0e-8        | Residual convergence tolerance. "property_map_D"        | int
   * | 0             | Material property index to use for diffusion coefficient
   * "property_map_q"        | int    | 1             | Material property index to
   * use for source. "property_map_sigma"    | int    | 2             | Material
   * property index to use for interaction coefficient.
   *
   * To set these options use the command chiSolverSetBasicOption() with
   * an option-name and value in the table above.
   *
   * ## Example
   * Example usage
   * \code
   * chiSolverSetBasicOption(phys1, "discretization_method", "PWLC")
   * \endcode
   */
  explicit Solver(const std::string& in_solver_name);
  virtual ~Solver();

  /**
   * Gets material properties various sources.
   */
  void GetMaterialProperties(const chi_mesh::Cell& cell,
                             int cell_dofs,
                             std::vector<double>& diffCoeff,
                             std::vector<double>& sourceQ,
                             std::vector<double>& sigmaa,
                             int group = 0,
                             int moment = 0);

  /**
   * Initialization of common to all solver types.
   */
  void InitializeCommonItems();

  void Initialize() override { Initialize(true); }

  /**
   * Initializes the diffusion solver using the PETSc library.
   */
  int Initialize(bool verbose);

  void Execute() override { ExecuteS(); }

  /**
   * Executes the diffusion solver using the PETSc library.
   */
  int ExecuteS(bool suppress_assembly = false, bool suppress_solve = false);

  /**
   * Assembles PWLC matrix for general cells.
   */
  void CFEM_Assemble_A_and_b(chi_mesh::Cell& cell, int group = 0);

  /**
   * Assembles PWLC matrix for polygon cells.
   */
  void PWLD_Assemble_A_and_b(const chi_mesh::Cell& cell, int component = 0);

  /**
   * Assembles PWLC matrix for polygon cells.
   */
  void PWLD_Assemble_b(const chi_mesh::Cell& cell, int component = 0);

  /**
   * Still searching for a reference for this.
   *
   * For Polygons:
   * Defined from paper  \n
   * Turcksin B, Ragusa J, "Discontinuous diffusion synthetic acceleration
   * for S_n transport on 2D arbitrary polygonal meshes", Journal of
   * Computational Physics 274, pg 356-369, 2014.\n
   * \n
   * Nv = Number of vertices. If Nv <= 4 then the perimeter parameter
   * should be replaced by edge length.
   */
  double HPerpendicular(const chi_mesh::Cell& cell,
                        const UnitIntegralContainer& fe_intgrl_values,
                        unsigned int f);

  /**
   * Given a global node index, returns the local cell-node it's associated on the referenced cell.
   */
  static uint64_t MapCellLocalNodeIDFromGlobalID(const chi_mesh::Cell& cell,
                                                 uint64_t node_global_id);

  /**
   * Given the face index on the current cell, finds the
   * corresponding face index on the adjacent cell.
   */
  static unsigned int
  MapCellFace(const chi_mesh::Cell& cur_cell, const chi_mesh::Cell& adj_cell, unsigned int f);

  /**
   * Update the field functions with the latest data.
   */
  void UpdateFieldFunctions();
};

} // namespace chi_diffusion
