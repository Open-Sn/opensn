#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "framework/math/Quadratures/angular_product_quadrature.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_group.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#define sc_int static_cast<int>

namespace lbs::common_lua_utils
{

/**Create a groupset.
\param SolverIndex int Handle to the solver for which the set is to be created.

##_

Example:
\code
gs0 = chiLBSCreateGroupset(phys1)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSCreateGroupset(lua_State* L)
{
  const std::string fname = "chiLBSCreateGroupset";
  // Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Create groupset
  lbs_solver.AddGroupset();

  lua_pushinteger(L, lbs_solver.Groupsets().back().id_);
  return 1;
}

/**Create a group.
\param SolverIndex int Handle to the solver for which the group
is to be created.
\param GroupId int Optional. If not supplied the next logical number
                             is assigned.

##_

Example:
\code
--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSCreateGroup(lua_State* L)
{
  const std::string fname = "chiLBSCreateGroup";
  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  // Get solver
  LuaCheckNumberValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);

  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Set Id
  int group_id = -1;
  if (num_args == 2)
  {
    LuaCheckNumberValue(fname, L, 2);
    group_id = lua_tointeger(L, 2);
  }

  // Create groupset
  lbs_solver.AddGroup(group_id);

  lua_pushinteger(L, lbs_solver.Groups().back().id_);
  return 1;
}

/**Adds a block of groups to a groupset.
\param SolverIndex int Handle to the solver for which the group
is to be created.
\param GroupsetIndex int Handle to the groupset to which the group is
 to be added.
\param FromIndex int From which group.
\param ToIndex int To which group.

##_

Example:
\code
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

chiLBSGroupsetAddGroups(phys1,cur_gs,0,15)
\endcode


\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetAddGroups(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetAddGroups";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 4) LuaPostArgAmountError(fname, 4, num_args);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int from = lua_tonumber(L, 3);
  int to = lua_tonumber(L, 4);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "chiLBSGroupsetAddGroups: Invalid handle to groupset\n";
    Chi::Exit(EXIT_FAILURE);
  }
  if (groupset == nullptr) throw std::runtime_error("chiLBSGroupsetAddGroups: Bad trouble.");

  // Add the groups
  if (to < from)
  {
    Chi::log.LogAllError() << "No groups added to groupset in chiLBSGroupsetAddGroups. "
                              "This is triggered when groups are added with the \"to\" "
                              "field being less than the \"from\" field.";
    Chi::Exit(EXIT_FAILURE);
  }

  for (unsigned k = from; k <= to; k++)
  {
    lbs::LBSGroup const* group = nullptr;
    // Check valid group
    try
    {
      group = &lbs_solver.Groups().at(k);
    }
    catch (const std::out_of_range& o)
    {
      Chi::log.LogAllError() << "chiLBSGroupsetAddGroups: Invalid group added to groupset\n";
      Chi::Exit(EXIT_FAILURE);
    }
    if (group == nullptr) throw std::runtime_error("chiLBSGroupsetAddGroups: Bad trouble.");

    groupset->groups_.push_back(*group);
  }
  return 0;
}

/**Sets the product quadrature used for the groupset
\param SolverIndex int Handle to the solver for which the group
is to be created.
\param GroupsetIndex int Handle to the groupset to which the group is
 to be added.
\param QuadratureIndex int Handle to the quadrature to be set for this
 groupset.


##_

Example:
\code
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)

chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetQuadrature(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetQuadrature";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int prquad_index = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset"
                              "in chiLBSGroupsetSetQuadrature.";
    Chi::Exit(EXIT_FAILURE);
  }
  if (groupset == nullptr) throw std::logic_error("chiLBSGroupsetSetQuadrature: Bad trouble");

  // Obtain pointer to quadrature
  std::shared_ptr<chi_math::AngularQuadrature> ang_quad;
  try
  {
    ang_quad = Chi::GetStackItemPtr(Chi::angular_quadrature_stack, prquad_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to Product Quadrature"
                              "in chiLBSGroupsetSetQuadrature. Handle provided: "
                           << prquad_index;
    Chi::Exit(EXIT_FAILURE);
  }

  groupset->quadrature_ = ang_quad;

  if (ang_quad->type_ == chi_math::AngularQuadratureType::ProductQuadrature)
  {
    auto prodquad = std::static_pointer_cast<chi_math::ProductQuadrature>(ang_quad);

    Chi::log.Log() << "Groupset " << grpset_index << " quadrature set to quadrature with "
                   << prodquad->azimu_ang_.size() << " azimuthal angles and "
                   << prodquad->polar_ang_.size() << " polar angles. ";
  }
  else if (ang_quad->type_ == chi_math::AngularQuadratureType::Arbitrary)
  {
    Chi::log.Log() << "Groupset " << grpset_index << " quadrature set to quadrature with "
                   << ang_quad->abscissae_.size() << " number of angles.";
  }
  else
    Chi::log.Log() << "Groupset " << grpset_index << " quadrature set unknown quadrature type";

  return 0;
}

/**Sets the the type of angle aggregation to use for this groupset.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param AggregationType int See AggregationType.

##_

###AggregationType
LBSGroupset.ANGLE_AGG_POLAR\n
 Use Polar angle aggregation. This is the default.\n\n

LBSGroupset.ANGLE_AGG_SINGLE\n
 Use Single angle aggregation.\n\n

LBSGroupset.ANGLE_AGG_AZIMUTHAL\n
 Use Azimuthal angle aggregation.\n\n

Example:
\code
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_POLAR)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetAngleAggregationType(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetAngleAggregationType";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int agg_type = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetAngleAggregationType";
    Chi::Exit(EXIT_FAILURE);
  }

  // Setting aggregation type
  if (agg_type == (int)lbs::AngleAggregationType::SINGLE)
    groupset->angleagg_method_ = lbs::AngleAggregationType::SINGLE;
  else if (agg_type == (int)lbs::AngleAggregationType::POLAR)
    groupset->angleagg_method_ = lbs::AngleAggregationType::POLAR;
  else if (agg_type == (int)lbs::AngleAggregationType::AZIMUTHAL)
    groupset->angleagg_method_ = lbs::AngleAggregationType::AZIMUTHAL;
  else
  {
    Chi::log.LogAllError() << "Invalid aggregation type to groupset " << grpset_index
                           << " in call to chiLBSGroupsetSetAngleAggregationType";
    Chi::Exit(EXIT_FAILURE);
  }

  Chi::log.Log() << "Groupset " << grpset_index << " Angle aggregation set to " << agg_type;

  return 0;
}

/**Sets the angle aggregation divisions
\param SolverIndex int Handle to the solver for which the group
is to be created.
\param GroupsetIndex int Handle to the groupset to which the group is
 to be added.
\param NumDiv int Number of divisions to use for the angle aggregation.

Note: by default polar aggregation will combine all polar angles in a hemisphere
 for a given azimuthal angleset. Therefore if there are 24 polar angles and
 4 azimuthal angles the default polar aggregation will create 8 anglesets
 (2 per quadrant to allow top and bottom hemisphere) and each angleset will have
the 12 polar angles associated with a hemisphere. When the number of divisions
is greater than 1 then the polar angles will be split into divisions. For
example if the number of divisions is 2 then more angleset will be created, this
time having 6 polar angles per angleset.

##_

Example:
\code
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetAngleAggDiv(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetAngleAggDiv";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int num_div = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetAngleAggDiv";
    Chi::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (num_div <= 0)
  {
    Chi::log.LogAllError() << "Invalid number of divisions "
                           << "in call to chiLBSGroupsetSetAngleAggDiv. Must be >= 1.";
  }

  groupset->master_num_ang_subsets_ = num_div;

  Chi::log.Log() << "Groupset " << grpset_index << " angle aggregation divisions "
                 << "set to " << num_div;

  return 0;
}

/**Sets the number of group-subsets to use for groupset. Default 1.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param NumDiv int Number of divisions of the groupset to use.

##_

Example:
\code
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetGroupSubsets(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetGroupSubsets";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int num_div = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetGroupSubsets";
    Chi::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (num_div <= 0)
  {
    Chi::log.LogAllError() << "Invalid number of subsets "
                           << "in call to chiLBSGroupsetSetGroupSubsets. Must be >= 1.";
  }

  groupset->master_num_grp_subsets_ = num_div;

  Chi::log.Log() << "Groupset " << grpset_index << " subset divisions "
                 << "set to " << num_div;

  return 0;
}

/**Sets the number of group-subsets to use for groupset. Default 1.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param IterativeMethod int Iteritve method identifier.

##_

### IterativeMethod
NPT_CLASSICRICHARDSON\n
Standard source iteration, without using PETSc.\n\n

NPT_CLASSICRICHARDSON_CYCLES\n
Standard source iteration, without using PETSc,
with cyclic dependency convergence.\n\n

NPT_GMRES\n
Legacy Generalized Minimal Residual formulation for iterations.\n\n

NPT_GMRES_CYCLES\n
Legacy Generalized Minimal Residual formulation for iterations with cyclic
dependency convergence.\n\n


KRYLOV_RICHARDSON\n
Richardson iteration.\n\n

KRYLOV_RICHARDSON_CYCLES\n
Richardson iteration with cyclic dependency convergence.\n\n

KRYLOV_GMRES\n
Generalized Minimal Residual method.\n\n

KRYLOV_GMRES_CYCLES\n
Generalized Minimal Residual method with cyclic dependency convergence.\n\n

KRYLOV_BICGSTAB\n
Biconjugate Gradient Stabilized method.\n\n

KRYLOV_BICGSTAB_CYCLES\n
Biconjugate Gradient Stabilized method with cyclic dependency convergence.\n\n


##_

### Notes on selecting iterative methods
The iterative methods NPT_CLASSICRICHARDSON, NPT_CLASSICRICHARDSON_CYCLES,
NPT_GMRES and NPT_GMRES_CYCLES are considered legacy. The NPT_GMRES and
NPT_GMRES_CYCLES are now considered deprecated with the inclusion of the
generalized Krylov iteration method-class (which supports all the options
prepended with KRYLOV_).

RICHARDSON is probably the least memory consuming but has the poorest
convergence rate.

GMRES generally has the best convergence rate but it builds a basis
comprising multiple solutions vectors, the amount of which is controlled via
the gmres-restart parameter, which can dramatically increase memory consumption.
GMRES restarts, i.e. the amount of iterations before the basis is destroyed and
restarted, influences both memory consumptions and convergence behavior, e.g.,
lower restart numbers generally lowers memory consumption but increases the
amount of iteration required to convergence.

The required memory and the computational time for one iteration with BiCGStab
is constant, i.e., the time and memory requirements do not increase with the
number of iterations as they do for restarted GMRES. BiCGStab uses approximately
the same amount of memory as GMRES uses for two iterations. Therefore, BiCGStab
typically uses less memory than GMRES. The convergence behavior of BiCGStab is
often more irregular than that of GMRES. Intermediate residuals can even be
orders of magnitude larger than the initial residual, which can affect the
numerical accuracy as well as the rate of convergence. If the algorithm detects
poor accuracy in the residual or the risk of stagnation, it restarts the
iterations with the current solution as the initial guess. In contrast to GMRES,
BiCGStab uses two matrix-vector multiplications each iteration (requiring two
transport sweeps). Also, when using the left-preconditioned BiCGStab, an
additional preconditioning step is required each iteration. That is,
left-preconditioned BiCGStab requires a total of three preconditioning steps in
each iteration. We generally apply Within-group Diffusion Synthetic Acceleration
(WGDSA) and Two-Grid Acceleration (TGDSA) as left-preconditioners and therefore
the total cost of these pre-conditioners will increase when using BiCGStab. Use
BiCGStab when you are running problem with a high scattering order (i.e., L is
large) because this will dramatically increase the GMRES basis.

Example:
\code
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_CLASSICRICHARDSON)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetIterativeMethod(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetIterativeMethod";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int iter_method = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetGroupSubsets";
    Chi::Exit(EXIT_FAILURE);
  }

  {
    using lbs::IterativeMethod;
    if (iter_method == sc_int(IterativeMethod::CLASSICRICHARDSON))
    {
      groupset->iterative_method_ = IterativeMethod::CLASSICRICHARDSON;
    }
    else if (iter_method == sc_int(IterativeMethod::CLASSICRICHARDSON_CYCLES))
    {
      groupset->allow_cycles_ = true;
      groupset->iterative_method_ = IterativeMethod::CLASSICRICHARDSON;
    }
    else if (iter_method == sc_int(IterativeMethod::GMRES))
    {
      throw std::invalid_argument(fname + "Deprecated iterative method GMRES, "
                                          "use KRYLOV_GMRES.");
    }
    else if (iter_method == sc_int(IterativeMethod::GMRES_CYCLES))
    {
      throw std::invalid_argument(fname + "Deprecated iterative method GMRES_CYCLES, "
                                          "use KRYLOV_GMRES_CYCLES.");
    }
    else if (iter_method == sc_int(IterativeMethod::KRYLOV_RICHARDSON))
    {
      groupset->iterative_method_ = IterativeMethod::KRYLOV_RICHARDSON;
    }
    else if (iter_method == sc_int(IterativeMethod::KRYLOV_RICHARDSON_CYCLES))
    {
      groupset->allow_cycles_ = true;
      groupset->iterative_method_ = IterativeMethod::KRYLOV_RICHARDSON;
    }
    else if (iter_method == sc_int(IterativeMethod::KRYLOV_GMRES))
    {
      groupset->iterative_method_ = IterativeMethod::KRYLOV_GMRES;
    }
    else if (iter_method == sc_int(IterativeMethod::KRYLOV_GMRES_CYCLES))
    {
      groupset->allow_cycles_ = true;
      groupset->iterative_method_ = IterativeMethod::KRYLOV_GMRES;
    }
    else if (iter_method == sc_int(IterativeMethod::KRYLOV_BICGSTAB))
    {
      groupset->iterative_method_ = IterativeMethod::KRYLOV_BICGSTAB;
    }
    else if (iter_method == sc_int(IterativeMethod::KRYLOV_BICGSTAB_CYCLES))
    {
      groupset->allow_cycles_ = true;
      groupset->iterative_method_ = IterativeMethod::KRYLOV_BICGSTAB;
    }
    else
    {
      Chi::log.LogAllError() << "Unsupported iterative method specified in call to "
                             << "chiLBSGroupsetSetIterativeMethod.";
      Chi::Exit(EXIT_FAILURE);
    }
  }

  return 0;
}

/**Sets the residual tolerance for the iterative method of the groupset.
 *
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param ResidualTol float residual tolerance (default 1.0e-6)

Note this tolerance also gets used for classic-richardson pointwise convergence
tolerance.

##_

Example:
\code
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetResidualTolerance(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetResidualTolerance";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  double resid_tol = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetGroupSubsets";
    Chi::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (resid_tol < 0)
  {
    Chi::log.LogAllError() << "Invalid residual tolerance specified. Must be greater >= 0.0";
    Chi::Exit(EXIT_FAILURE);
  }

  groupset->residual_tolerance_ = resid_tol;

  char buff[100];
  snprintf(buff, 100, "%.4e", resid_tol);

  Chi::log.Log() << "Groupset " << grpset_index << " residual tolerance "
                 << "set to " << buff;

  return 0;
}

/**Sets the maximum number of iterations for the groupset iterative method.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param Numiter int Maximum number of iterations. Default 1000.

##_

Example:
\code
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,200)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetMaxIterations(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetMaxIterations";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int num_iter = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetMaxIterations";
    Chi::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (num_iter < 0)
  {
    Chi::log.LogAllError() << "Invalid number of iterations "
                           << "in call to chiLBSGroupsetSetMaxIterations. Must be >= 0.";
  }

  groupset->max_iterations_ = num_iter;

  Chi::log.Log() << "Groupset " << grpset_index << " max # iterations "
                 << "set to " << num_iter;

  return 0;
}

/**Sets the restart interval for GMRES if applied to the groupset.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param Intvl int Interval to use for GMRES restarts. Default 30.

##_

Example:
\code
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,15)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetGMRESRestartIntvl(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetGMRESRestartIntvl";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int restart_intvl = lua_tonumber(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetGMRESRestartIntvl";
    Chi::Exit(EXIT_FAILURE);
  }

  // Bounds checking
  if (restart_intvl < 2)
  {
    Chi::log.LogAllError() << "Invalid GMRES restart interval specified "
                           << "in call to chiLBSGroupsetSetGMRESRestartIntvl. Must be >= 3.";
  }

  groupset->gmres_restart_intvl_ = restart_intvl;

  Chi::log.Log() << "Groupset " << grpset_index << " GMRES restart interval set to "
                 << "set to " << restart_intvl;

  return 0;
}

/**Enables or disables the printing of a sweep log.
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param flag bool Flag indicating whether to print sweep log. Default false.

##_

Example:
\code
chiLBSGroupsetSetEnableSweepLog(phys1,cur_gs,true)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetEnableSweepLog(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetEnableSweepLog";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  bool log_flag = lua_toboolean(L, 3);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetEnableSweepLog";
    Chi::Exit(EXIT_FAILURE);
  }

  groupset->log_sweep_events_ = log_flag;

  Chi::log.Log() << "Groupset " << grpset_index << " flag for writing sweep log "
                 << "set to " << log_flag;

  return 0;
}

/**Sets the Within-Group Diffusion Synthetic Acceleration parameters
 * for this groupset. If this call is being made then it is assumed
 * WGDSA is being applied.
 *
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param MaxIters int Maximum amount of iterations to use for WGDSA solvers.
                    Default 30.
\param ResTol float Residual tolerance to use for the WGDSA solve.

\param Verbose bool Optional flag indicating verbose output of WGDSA.
                    Default false.
\param PETSCString char Optional. Options string to be inserted
                        during initialization.



##_

Example:
\code
petsc_options =                  " -pc_hypre_boomeramg_strong_threshold 0.8"
petsc_options = petsc_options .. " -pc_hypre_boomeramg_max_levels 25"
chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false,petsc_options)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetWGDSA(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetWGDSA";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args < 4) LuaPostArgAmountError(fname, 4, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  LuaCheckNilValue(fname, L, 4);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int max_iters = lua_tonumber(L, 3);
  double resid_tol = lua_tonumber(L, 4);
  bool verbose = false;
  const char* petsc_string = "";

  if (num_args >= 5) verbose = lua_toboolean(L, 5);

  if (num_args == 6) petsc_string = lua_tostring(L, 6);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetWGDSA";
    Chi::Exit(EXIT_FAILURE);
  }

  groupset->apply_wgdsa_ = true;
  groupset->wgdsa_max_iters_ = max_iters;
  groupset->wgdsa_tol_ = resid_tol;
  groupset->wgdsa_verbose_ = verbose;
  groupset->wgdsa_string_ = std::string(petsc_string);

  Chi::log.Log() << "Groupset " << grpset_index << " set to apply WGDSA with " << max_iters
                 << " maximum iterations and a tolerance of " << resid_tol
                 << ". PETSc-string: " << std::string(petsc_string);

  return 0;
}

/**Sets the Two-Grid Diffusion Synthetic Acceleration parameters
 * for this groupset. If this call is being made then it is assumed
 * TGDSA is being applied.
 *
\param SolverIndex int Handle to the solver for which the group
is to be created.

\param GroupsetIndex int Index to the groupset to which this function should
                         apply
\param MaxIters int Maximum amount of iterations to use for TGDSA solvers.
                    Default 30.
\param ResTol float Residual tolerance to use for the TGDSA solve.

\param Verbose bool Optional flag indicating verbose output of TGDSA.
                    Default false.
\param PETSCString char Optional. Options string to be inserted
                        during initialization.



##_

Example:
\code
petsc_options =                  " -pc_hypre_boomeramg_strong_threshold 0.8"
petsc_options = petsc_options .. " -pc_hypre_boomeramg_max_levels 25"
chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false,petsc_options)
\endcode

\ingroup LuaLBSGroupsets
*/
int
chiLBSGroupsetSetTGDSA(lua_State* L)
{
  const std::string fname = "chiLBSGroupsetSetTGDSA";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args < 4) LuaPostArgAmountError(fname, 4, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  LuaCheckNilValue(fname, L, 4);
  int solver_handle = lua_tonumber(L, 1);
  int grpset_index = lua_tonumber(L, 2);
  int max_iters = lua_tonumber(L, 3);
  double resid_tol = lua_tonumber(L, 4);
  bool verbose = false;
  const char* petsc_string = "";

  if (num_args >= 5) verbose = lua_toboolean(L, 5);

  if (num_args == 6) petsc_string = lua_tostring(L, 6);

  // Get pointer to solver
  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    Chi::log.LogAllError() << "Invalid handle to groupset "
                           << "in call to chiLBSGroupsetSetTGDSA";
    Chi::Exit(EXIT_FAILURE);
  }

  groupset->apply_tgdsa_ = true;
  groupset->tgdsa_max_iters_ = max_iters;
  groupset->tgdsa_tol_ = resid_tol;
  groupset->tgdsa_verbose_ = verbose;
  groupset->tgdsa_string_ = std::string(petsc_string);

  Chi::log.Log() << "Groupset " << grpset_index << " set to apply TGDSA with " << max_iters
                 << " maximum iterations and a tolerance of " << resid_tol
                 << ". PETSc-string: " << std::string(petsc_string);

  return 0;
}

} // namespace lbs::common_lua_utils
