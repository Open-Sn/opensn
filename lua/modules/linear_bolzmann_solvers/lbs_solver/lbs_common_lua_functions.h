#pragma once

#include "framework/lua.h"

namespace opensnlua::lbs
{
/**Set LBS property.
 * \param SolverIndex int Handle to the solver for which the set is to be created.
 * \param PropertyIndex int Code for a specific property.
 *
 *
 * ##_
 *
 * ###PropertyIndex
 * DISCRETIZATION_METHOD\n
 *  Discretization method.\n\n
 *
 * BOUNDARY_CONDITION\n
 *  Boundary condition type. See BoundaryIdentify.\n\n
 *
 * SCATTERING_ORDER\n
 *  Defines the level of harmonic expansion for the scattering source.Default 1.
 *  Expects to be followed by an integer.\n\n
 *
 * SWEEP_EAGER_LIMIT\n
 *  The eager limit to be used in message size during sweep initialization.
 *  This expects to be followed by a size in bytes (Max 64,0000).Default 32,000.
 *  See note below.\n\n
 *
 * READ_RESTART_DATA\n
 *  Indicates the reading of restart data from restart file.
 *  The value can be followed by two
 *  optional strings. The first is the folder name which can be relative or
 *  absolute, and the second is the file base name. These are defaulted to
 *  "YRestart" and "restart" respectively.\n\n
 *
 * SAVE_ANGULAR_FLUX\n
 * Sets the flag for saving the angular flux. Expects to be followed by true/false.
 * [Default=false]\n\n
 *
 * USE_SOURCE_MOMENTS\n
 *  Flag for using a vector of source moments instead the regular material/boundary
 *   source. Default false. This expects
 *  to be followed by a boolean.\n\n
 *
 * VERBOSE_INNER_ITERATIONS\n
 *  Flag for printing inner iteration information. This is primarily used
 *  for printing information related to group-set-level iterative methods.
 *  Default true. Expects to be followed by a boolean.\n\n
 *
 * VERBOSE_OUTER_ITERATIONS\n
 *  Flag for printing outer iteration information. This is primarily used
 *  for printing information aggregated over group sets such as k-eigenvalue
 *  iterations. Default true. Expects to be followed by a boolean.\n\n
 *
 * USE_PRECURSORS\n
 *  Flag for using delayed neutron precursors. Default false. This expects
 *  to be followed by a boolean.\n\n
 *
 * \code
 * LBSSetProperty(phys1,READ_RESTART_DATA,"YRestart1")
 * \endcode
 *
 * WRITE_RESTART_DATA\n
 *  Indicates the writing of restart data to restart files.
 *  The value can be followed by two optional strings and a number
 *  optional strings. The first string is the folder name which can be relative or
 *  absolute, and the second string is the file base name. The number is the time
 *  interval (in minutes) for a restart write to be triggered (apart from GMRES
 *  restarts and the conclusion of groupset completions) .These are defaulted to
 *  "YRestart", "restart" and 30 minutes respectively.\n\n
 *
 * \code
 * LBSSetProperty(phys1,WRITE_RESTART_DATA,"YRestart1","restart",1)
 * \endcode
 *
 * ###Discretization methods
 *  PWLD = Piecewise Linear Finite Element.\n
 *
 * ###BoundaryIdentify
 * This value follows the argument BOUNDARY_CONDITION and identifies which
 * boundary is under consideration. Right now only boundaries aligned with
 * cartesian axes are considered. Followed by LBSBoundaryType.\n
 * XMAX = Right boundary \n
 * XMIN = Left boundary \n
 * YMAX = Front boundary \n
 * YMIN = Back boundary \n
 * ZMAX = Top boundary \n
 * ZMIN = Bottom boundary \n
 *
 * ###LBSBoundaryType
 * Specifies the type of boundary. Depending on the type this argument needs
 * to be followed by one or more values. Note: By default all boundaries are
 * type VACUUM.\n
 * \n
 * LBSBoundaryTypes.VACUUM\n
 * Specifies a vaccuum boundary condition. It is not followed by any value.\n
 * \code
 * LBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,
 *                       LBSBoundaryTypes.VACUUM);
 * \endcode
 * \n
 * LBSBoundaryTypes.INCIDENT_ISOTROPIC\n
 * Incident isotropic flux. This argument needs to be followed by a lua table
 * index 1 to G where G is the amount of energy groups. Note internally this
 * is mapped as 0 to G-1.\n
 * \code
 * bsrc={}
 * for g=1,num_groups do
 *     bsrc[g] = 0.0
 * end
 * bsrc[1] = 1.0
 * LBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,
 *                       LBSBoundaryTypes.INCIDENT_ISOTROPIC, bsrc);
 * \endcode
 * \n
 * LBSBoundaryTypes.REFLECTING\n
 * Reflecting boundary condition. Beware, when opposing reflecting boundary
 * conditions are used this enduces a cyclic dependency which will increase the
 * iteration convergence behavior.\n
 * \code
 * LBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,
 *                       LBSBoundaryTypes.REFLECTING);
 * \endcode
 * \n
 * LBSBoundaryTypes.INCIDENT_ANISTROPIC_HETEROGENEOUS\n
 * Expects to be followed by the name of a lua function. The lua function will get
 * called with the following parameters:
 * ```
 * size_t        cell_global_id,
 * int           cell_material_id,
 * unsigned int  face_index,
 * unsigned int  face_node_index,
 * const Vector3& face_node_location,
 * const Vector3& face_node_normal,
 * const std::vector<int>& quadrature_angle_indices,
 * const std::vector<Vector3>& quadrature_angle_vectors,
 * const std::vector<std::pair<double,double>>& quadrature_phi_theta_angles,
 * const std::vector<int>& group_indices,
 * double evaluation_time;
 * ```
 * and must return a 1D array of data-values ordered first by angle index, then
 * by group index, e.g., n0g0, n0g1, n0g2, n1g0, n1g1, n1g2, etc.
 *
 * Example lua function:
 * \code
 * function luaBoundaryFunctionA(cell_global_id,
 *                               material_id,
 *                               location,
 *                               normal,
 *                               quadrature_angle_indices,
 *                               quadrature_angle_vectors,
 *                               quadrature_phi_theta_angles,
 *                               group_indices)
 *     num_angles = rawlen(quadrature_angle_vectors)
 *     num_groups = rawlen(group_indices)
 *     psi = {}
 *     dof_count = 0
 *
 *     for ni=1,num_angles do
 *         omega = quadrature_angle_vectors[ni]
 *         phi_theta = quadrature_phi_theta_angles[ni]
 *         for gi=1,num_groups do
 *             g = group_indices[gi]
 *
 *             value = 1.0
 *
 *             dof_count = dof_count + 1
 *             psi[dof_count] = value
 *         end
 *     end
 *
 *     return psi
 * end
 *
 * LBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,
 *                       LBSBoundaryTypes.INCIDENT_ANISTROPIC_HETEROGENEOUS,
 *                       "luaBoundaryFunctionA");
 * \endcode
 *
 *
 * ###Note on the Eager limit
 * The eager limit is the message size limit before which non-blocking MPI send
 * calls will execute without waiting for a matching receive call. The limit is
 * platform dependent but in general 64 kb. Some systems have 32 kb as a limit
 * and therefore we use that as a default limit in OpenSn. There is a fine
 * interplay between message size and the shear amount of messages that will be
 * sent. In general smaller messages tend to be more efficient, however, when
 * there are too many small messages being sent around the communication system
 * on the given platform will start to suffer. One can gain a small amount of
 * parallel efficiency by lowering this limit, however, there is a point where
 * the parallel efficiency will actually get worse so use with caution.
 *
 * \ingroup LBSLuaFunctions
 */
int LBSSetProperty(lua_State* L);

/**Create a groupset.
 * \param SolverIndex int Handle to the solver for which the set is to be created.
 *
 * ##_
 *
 * Example:
 * \code
 * gs0 = LBSCreateGroupset(phys1)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSCreateGroupset(lua_State* L);

/**Create a group.
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 * \param GroupId int Optional. If not supplied the next logical number
 *                              is assigned.
 *
 * ##_
 *
 * Example:
 * \code
 * --========== Groups
 * grp = {}
 * for g=1,num_groups do
 *     grp[g] = LBSCreateGroup(phys1)
 * end
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSCreateGroup(lua_State* L);

/**Adds a block of groups to a groupset.
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 * \param GroupsetIndex int Handle to the groupset to which the group is
 *  to be added.
 * \param FromIndex int From which group.
 * \param ToIndex int To which group.
 *
 * ##_
 *
 * Example:
 * \code
 * grp = {}
 * for g=1,num_groups do
 *     grp[g] = LBSCreateGroup(phys1)
 * end
 *
 * LBSGroupsetAddGroups(phys1,cur_gs,0,15)
 * \endcode
 *
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetAddGroups(lua_State* L);

/**Sets the product quadrature used for the groupset
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 * \param GroupsetIndex int Handle to the groupset to which the group is
 *  to be added.
 * \param QuadratureIndex int Handle to the quadrature to be set for this
 *  groupset.
 *
 *
 * ##_
 *
 * Example:
 * \code
 * pquad0 = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
 *
 * LBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetQuadrature(lua_State* L);

/**Sets the the type of angle aggregation to use for this groupset.
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param AggregationType int See AggregationType.
 *
 * ##_
 *
 * ###AggregationType
 * LBSGroupset.ANGLE_AGG_POLAR\n
 *  Use Polar angle aggregation. This is the default.\n\n
 *
 * LBSGroupset.ANGLE_AGG_SINGLE\n
 *  Use Single angle aggregation.\n\n
 *
 * LBSGroupset.ANGLE_AGG_AZIMUTHAL\n
 *  Use Azimuthal angle aggregation.\n\n
 *
 * Example:
 * \code
 * LBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_POLAR)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetAngleAggregationType(lua_State* L);

/**Sets the angle aggregation divisions
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 * \param GroupsetIndex int Handle to the groupset to which the group is
 *  to be added.
 * \param NumDiv int Number of divisions to use for the angle aggregation.
 *
 * Note: by default polar aggregation will combine all polar angles in a hemisphere
 *  for a given azimuthal angleset. Therefore if there are 24 polar angles and
 *  4 azimuthal angles the default polar aggregation will create 8 anglesets
 *  (2 per quadrant to allow top and bottom hemisphere) and each angleset will have
 * the 12 polar angles associated with a hemisphere. When the number of divisions
 * is greater than 1 then the polar angles will be split into divisions. For
 * example if the number of divisions is 2 then more angleset will be created, this
 * time having 6 polar angles per angleset.
 *
 * ##_
 *
 * Example:
 * \code
 * LBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetAngleAggDiv(lua_State* L);

/**Sets the number of group-subsets to use for groupset. Default 1.
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param NumDiv int Number of divisions of the groupset to use.
 *
 * ##_
 *
 * Example:
 * \code
 * LBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetGroupSubsets(lua_State* L);

/**Sets the number of group-subsets to use for groupset. Default 1.
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param IterativeMethod int Iteritve method identifier.
 *
 * ##_
 *
 * ### IterativeMethod
 * NPT_CLASSICRICHARDSON\n
 * Standard source iteration, without using PETSc.\n\n
 *
 * NPT_CLASSICRICHARDSON_CYCLES\n
 * Standard source iteration, without using PETSc,
 * with cyclic dependency convergence.\n\n
 *
 * NPT_GMRES\n
 * Legacy Generalized Minimal Residual formulation for iterations.\n\n
 *
 * NPT_GMRES_CYCLES\n
 * Legacy Generalized Minimal Residual formulation for iterations with cyclic
 * dependency convergence.\n\n
 *
 *
 * KRYLOV_RICHARDSON\n
 * Richardson iteration.\n\n
 *
 * KRYLOV_RICHARDSON_CYCLES\n
 * Richardson iteration with cyclic dependency convergence.\n\n
 *
 * KRYLOV_GMRES\n
 * Generalized Minimal Residual method.\n\n
 *
 * KRYLOV_GMRES_CYCLES\n
 * Generalized Minimal Residual method with cyclic dependency convergence.\n\n
 *
 * KRYLOV_BICGSTAB\n
 * Biconjugate Gradient Stabilized method.\n\n
 *
 * KRYLOV_BICGSTAB_CYCLES\n
 * Biconjugate Gradient Stabilized method with cyclic dependency convergence.\n\n
 *
 *
 * ##_
 *
 * ### Notes on selecting iterative methods
 * The iterative methods NPT_CLASSICRICHARDSON, NPT_CLASSICRICHARDSON_CYCLES,
 * NPT_GMRES and NPT_GMRES_CYCLES are considered legacy. The NPT_GMRES and
 * NPT_GMRES_CYCLES are now considered deprecated with the inclusion of the
 * generalized Krylov iteration method-class (which supports all the options
 * prepended with KRYLOV_).
 *
 * RICHARDSON is probably the least memory consuming but has the poorest
 * convergence rate.
 *
 * GMRES generally has the best convergence rate but it builds a basis
 * comprising multiple solutions vectors, the amount of which is controlled via
 * the gmres-restart parameter, which can dramatically increase memory consumption.
 * GMRES restarts, i.e. the amount of iterations before the basis is destroyed and
 * restarted, influences both memory consumptions and convergence behavior, e.g.,
 * lower restart numbers generally lowers memory consumption but increases the
 * amount of iteration required to convergence.
 *
 * The required memory and the computational time for one iteration with BiCGStab
 * is constant, i.e., the time and memory requirements do not increase with the
 * number of iterations as they do for restarted GMRES. BiCGStab uses approximately
 * the same amount of memory as GMRES uses for two iterations. Therefore, BiCGStab
 * typically uses less memory than GMRES. The convergence behavior of BiCGStab is
 * often more irregular than that of GMRES. Intermediate residuals can even be
 * orders of magnitude larger than the initial residual, which can affect the
 * numerical accuracy as well as the rate of convergence. If the algorithm detects
 * poor accuracy in the residual or the risk of stagnation, it restarts the
 * iterations with the current solution as the initial guess. In contrast to GMRES,
 * BiCGStab uses two matrix-vector multiplications each iteration (requiring two
 * transport sweeps). Also, when using the left-preconditioned BiCGStab, an
 * additional preconditioning step is required each iteration. That is,
 * left-preconditioned BiCGStab requires a total of three preconditioning steps in
 * each iteration. We generally apply Within-group Diffusion Synthetic Acceleration
 * (WGDSA) and Two-Grid Acceleration (TGDSA) as left-preconditioners and therefore
 * the total cost of these pre-conditioners will increase when using BiCGStab. Use
 * BiCGStab when you are running problem with a high scattering order (i.e., L is
 * large) because this will dramatically increase the GMRES basis.
 *
 * Example:
 * \code
 * LBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_CLASSICRICHARDSON)
 * LBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetIterativeMethod(lua_State* L);

/**Sets the residual tolerance for the iterative method of the groupset.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param ResidualTol float residual tolerance (default 1.0e-6)
 *
 * Note this tolerance also gets used for classic-richardson pointwise convergence
 * tolerance.
 *
 * ##_
 *
 * Example:
 * \code
 * LBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetResidualTolerance(lua_State* L);

/**Sets the maximum number of iterations for the groupset iterative method.
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param Numiter int Maximum number of iterations. Default 1000.
 *
 * ##_
 *
 * Example:
 * \code
 * LBSGroupsetSetMaxIterations(phys1,cur_gs,200)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetMaxIterations(lua_State* L);

/**Sets the restart interval for GMRES if applied to the groupset.
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param Intvl int Interval to use for GMRES restarts. Default 30.
 *
 * ##_
 *
 * Example:
 * \code
 * LBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,15)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetGMRESRestartIntvl(lua_State* L);

/**Enables or disables the printing of a sweep log.
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param flag bool Flag indicating whether to print sweep log. Default false.
 *
 * ##_
 *
 * Example:
 * \code
 * LBSGroupsetSetEnableSweepLog(phys1,cur_gs,true)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetEnableSweepLog(lua_State* L);

/**Sets the Within-Group Diffusion Synthetic Acceleration parameters
 * for this groupset. If this call is being made then it is assumed
 * WGDSA is being applied.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param MaxIters int Maximum amount of iterations to use for WGDSA solvers.
 *                     Default 30.
 * \param ResTol float Residual tolerance to use for the WGDSA solve.
 *
 * \param Verbose bool Optional flag indicating verbose output of WGDSA.
 *                     Default false.
 * \param PETSCString char Optional. Options string to be inserted
 *                         during initialization.
 *
 *
 *
 * ##_
 *
 * Example:
 * \code
 * petsc_options =                  " -pc_hypre_boomeramg_strong_threshold 0.8"
 * petsc_options = petsc_options .. " -pc_hypre_boomeramg_max_levels 25"
 * LBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false,petsc_options)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetWGDSA(lua_State* L);

/**Sets the Two-Grid Diffusion Synthetic Acceleration parameters
 * for this groupset. If this call is being made then it is assumed
 * TGDSA is being applied.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 * \param MaxIters int Maximum amount of iterations to use for TGDSA solvers.
 *                     Default 30.
 * \param ResTol float Residual tolerance to use for the TGDSA solve.
 *
 * \param Verbose bool Optional flag indicating verbose output of TGDSA.
 *                     Default false.
 * \param PETSCString char Optional. Options string to be inserted
 *                         during initialization.
 *
 *
 *
 * ##_
 *
 * Example:
 * \code
 * petsc_options =                  " -pc_hypre_boomeramg_strong_threshold 0.8"
 * petsc_options = petsc_options .. " -pc_hypre_boomeramg_max_levels 25"
 * LBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false,petsc_options)
 * \endcode
 *
 * \ingroup LuaLBSGroupsets
 */
int LBSGroupsetSetTGDSA(lua_State* L);

/**Obtains a list of field functions, related only to scalar flux,
 * from the transport solver.
 *
 * \param SolverIndex int Handle to the solver for which the list is to be
 * obtained.
 *
 * \return Pair Table and count. Returns an array of handles and the amount of
 * elements in it (indexed from 1). \ingroup LBSLuaFunctions \author Jan
 */
int LBSGetScalarFieldFunctionList(lua_State* L);

/**Writes the angular fluxes of a LBS groupset to file.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 *
 */
int chiLBSWriteGroupsetAngularFlux(lua_State* L);

/**Reads the angular fluxes of a LBS groupset from a file.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 *
 */
int chiLBSReadGroupsetAngularFlux(lua_State* L);

/**Writes the flux-moments of a LBS solution to file (phi_old_local).
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 *
 */
int LBSWriteFluxMoments(lua_State* L);

/**Creates scattered source-moments, based on a LBS solution, and writes them
 * to file.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 *
 */
int chiLBSCreateAndWriteSourceMoments(lua_State* L);

/**Reads flux-moments from a file and creates a scattering source from these
 * moments to be used instead of a regular material/boundary source.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 *
 * \param single_file_flag bool (Optional) Flag indicating that the file is a
 *                              single stand-alone file. The file_base will then
 *                              be used without adding the location-id, but still
 *                              with the ".data" appended. Default: false.
 *
 */
int chiLBSReadFluxMomentsAndMakeSourceMoments(lua_State* L);

/**Reads the source-moments from a file to a specific
 * ext_src_moments_local-vector
 *  * to be used instead of a regular material/boundary source.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 *
 * \param single_file_flag bool (Optional) Flag indicating that the file is a
 *                              single stand-alone file. The file_base will then
 *                              be used without adding the location-id, but still
 *                              with the ".data" appended. Default: false.
 */
int chiLBSReadSourceMoments(lua_State* L);

/**Reads flux-moments from a file to phi_old_local (the initial flux solution).
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 *
 * \param single_file_flag bool (Optional) Flag indicating that the file is a
 *                              single stand-alone file. The file_base will then
 *                              be used without adding the location-id, but still
 *                              with the ".data" appended. Default: false.
 */
int chiLBSReadFluxMoments(lua_State* L);

/**Computes and returns the fission rate.
 *
 * \param SolverIndex int Handle to the solver maintaining the information.
 * \param OldNewOption string "NEW" or "OLD". For steady state solvers, the
 *                            "OLD" option would give the fission rate for
 *                            the previous iterate. [Default="NEW"]
 *
 * \return double The fission rate.
 *
 * \ingroup LBSLuaFunctions
 * \author Jan
 */
int chiLBSComputeFissionRate(lua_State* L);

/**Initializes or reinitializes the materials. This normally happens
 * automatically during solver initialization but if the user wants to
 * swap/change XSs during the run then this will allow the material structures
 * to now deal with the new/changed materials.
 *
 * \param SolverIndex int Handle to the solver maintaining the information.
 *
 * \ingroup LBSLuaFunctions
 * \author Jan
 */
int chiLBSInitializeMaterials(lua_State* L);

/**Adds a point source to an LBS solver.
 * \param SolverIndex int Handle to the solver.
 * \param Location_x double X-location.
 * \param Location_y double Y-location.
 * \param Location_z double Z-location.
 * \param Strength table Source strength as a multigroup vector.
 *
 *  \ingroup LBSLuaFunctions
 */
int LBSAddPointSource(lua_State* L);

/**Clears all the point sources from the solver. This is mostly
 * useful for adjoint response calculations.
 * \param SolverIndex int Handle to the solver.
 *
 *  \ingroup LBSLuaFunctions
 */
int chiLBSClearPointSources(lua_State* L);

} // namespace opensnlua::lbs
