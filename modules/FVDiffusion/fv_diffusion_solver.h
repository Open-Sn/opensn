#pragma once

#include "opensn/framework/physics/SolverBase/chi_solver.h"
#include "opensn/framework/math/PETScUtils/petsc_utils.h"

#include "opensn/modules/FVDiffusion/fv_diffusion_bndry.h"
#include "opensn/framework/utils/chi_timer.h"

#include "opensn/framework/console/chi_console.h"

#include "opensn/framework/mesh/chi_mesh.h"

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

  Vec x_ = nullptr; // approx solution
  Vec b_ = nullptr; // RHS
  Mat A_ = nullptr; // linear system matrix

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
  static double
  CallLua_iXYZFunction(lua_State* L, const std::string&, int, const chi_mesh::Vector3&);
#endif

  void UpdateFieldFunctions();
};

} // namespace fv_diffusion
