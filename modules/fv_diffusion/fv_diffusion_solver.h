#pragma once

#include "framework/physics/solver_base/solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"

#include "modules/fv_diffusion/fv_diffusion_bndry.h"
#include "framework/utils/timer.h"

#include "framework/console/console.h"

#include "framework/mesh/mesh.h"

#include <map>

// forward declaration
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

namespace fv_diffusion
{
/** FV diffusion solver
 *
 */
class Solver : public chi_physics::Solver
{
public:
  chi_mesh::MeshContinuumPtr grid_ptr_ = nullptr;

  chi_math::SDMPtr sdm_ptr_ = nullptr;

  size_t num_local_dofs_ = 0;
  size_t num_globl_dofs_ = 0;

  /// approx solution
  Vec x_ = nullptr;
  /// RHS
  Vec b_ = nullptr;
  /// linear system matrix
  Mat A_ = nullptr;

  typedef std::pair<fv_diffusion::BoundaryType, std::vector<double>> BoundaryInfo;
  typedef std::map<std::string, BoundaryInfo> BoundaryPreferences;
  BoundaryPreferences boundary_preferences_;
  std::map<uint64_t, Boundary> boundaries_;

  explicit Solver(const std::string& in_solver_name);
  ~Solver() override;

  // void Initialize() override;
  void Initialize() override;
  void Execute() override;

#ifdef OPENSN_WITH_LUA
  /**Calls a lua function with xyz coordinates.
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

  /**Updates the field functions with the latest data.*/
  void UpdateFieldFunctions();
};

} // namespace fv_diffusion
