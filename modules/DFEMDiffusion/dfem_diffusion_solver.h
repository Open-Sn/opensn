#pragma once

#include "framework/physics/SolverBase/chi_solver.h"
#include "framework/math/PETScUtils/petsc_utils.h"
#include "modules/DFEMDiffusion/dfem_diffusion_bndry.h"
#include "framework/utils/chi_timer.h"
#include "framework/console/chi_console.h"
#include "framework/math/UnknownManager/unknown_manager.h"
#include "framework/mesh/chi_mesh.h"
#include <map>

namespace chi_mesh
{
class MeshContinuum;
typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
} // namespace chi_mesh

namespace chi_math
{
class SpatialDiscretization;
typedef std::shared_ptr<SpatialDiscretization> SDMPtr;
} // namespace chi_math

namespace dfem_diffusion
{

/**
 * DFEM diffusion solver
 */
class Solver : public chi_physics::Solver
{
public:
  chi_mesh::MeshContinuumPtr grid_ptr_ = nullptr;

  chi_math::SDMPtr sdm_ptr_ = nullptr;

  size_t num_local_dofs_ = 0;
  size_t num_globl_dofs_ = 0;

  std::vector<double> field_;

  /// approx solution
  Vec x_ = nullptr;
  /// RHS
  Vec b_ = nullptr;
  /// linear system matrix
  Mat A_ = nullptr;

  typedef std::pair<BoundaryType, std::vector<double>> BoundaryInfo;
  typedef std::map<std::string, BoundaryInfo> BoundaryPreferences;
  BoundaryPreferences boundary_preferences_;
  std::map<uint64_t, Boundary> boundaries_;

  explicit Solver(const std::string& in_solver_name);
  ~Solver() override;

  void Initialize() override;
  void Execute() override;

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
  double HPerpendicular(const chi_mesh::Cell& cell, unsigned int f);

  /**
   * Maps a face, in a discontinuous sense, using the spatial discretization.
   */
  int MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                      const chi_mesh::Cell& adj_cell,
                      const std::vector<chi_mesh::Vector3>& cc_node_locs,
                      const std::vector<chi_mesh::Vector3>& ac_node_locs,
                      size_t ccf,
                      size_t acf,
                      size_t ccfi,
                      double epsilon = 1.0e-12);

#ifdef OPENSN_WITH_LUA
  /**
   * Calls a lua function with xyz coordinates.
   * \param L The lua state.
   * \param lua_func_name The name used to define this lua function in the lua
   *                      state.
   * \param imat The material ID of the cell
   * \param xyz The xyz coordinates of the point where the function is called.
   *
   * \return The function evaluation.*/
  static double
  CallLua_iXYZFunction(lua_State* L, const std::string&, int, const chi_mesh::Vector3&);
#endif

  /**
   * Updates the field functions with the latest data.
   */
  void UpdateFieldFunctions();
};

} // namespace dfem_diffusion
