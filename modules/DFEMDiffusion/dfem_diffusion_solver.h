#ifndef DFEM_DIFFUSION_SOLVER_H
#define DFEM_DIFFUSION_SOLVER_H

#include "physics/SolverBase/chi_solver.h"
#include "math/PETScUtils/petsc_utils.h"

#include "dfem_diffusion_bndry.h"
#include "utils/chi_timer.h"

#include "console/chi_console.h"
#include "math/UnknownManager/unknown_manager.h"

#include "mesh/chi_mesh.h"

#include <map>

// forward declaration
namespace chi_mesh
{
class MeshContinuum; 
typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
}
namespace chi_math
{
class SpatialDiscretization; 
typedef std::shared_ptr<SpatialDiscretization> SDMPtr ;
}

namespace dfem_diffusion
{
/** DFEM diffusion solver
 * 
*/
class Solver : public chi_physics::Solver
{
public:
  chi_mesh::MeshContinuumPtr grid_ptr_ = nullptr;

  chi_math::SDMPtr sdm_ptr_ = nullptr;

  size_t num_local_dofs_ = 0;
  size_t num_globl_dofs_ = 0;

  std::vector<double> field_;

  Vec            x_ = nullptr;            // approx solution
  Vec            b_ = nullptr;            // RHS
  Mat            A_ = nullptr;            // linear system matrix

  typedef std::pair<BoundaryType,std::vector<double>> BoundaryInfo;
  typedef std::map<std::string, BoundaryInfo> BoundaryPreferences;
  BoundaryPreferences      boundary_preferences_;
  std::map<uint64_t, Boundary>    boundaries_;

  explicit Solver(const std::string& in_solver_name);
  ~Solver() override;

  // void Initialize() override;
  void Initialize() override;
  void Execute() override;

  double HPerpendicular(const chi_mesh::Cell& cell, unsigned int f);

  int MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                      const chi_mesh::Cell& adj_cell,
                      const std::vector<chi_mesh::Vector3>& cc_node_locs,
                      const std::vector<chi_mesh::Vector3>& ac_node_locs,
                      size_t ccf, size_t acf,
                      size_t ccfi,
                      double epsilon=1.0e-12);

  static double CallLua_iXYZFunction(lua_State* L,
                                     const std::string&,
                                     int,
                                     const chi_mesh::Vector3&);

  void UpdateFieldFunctions();
};

} // namespace dfem_diffusion


#endif //DFEM_DIFFUSION_SOLVER_H

