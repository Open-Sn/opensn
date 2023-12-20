#pragma once

namespace opensn
{
class InputParameters;
}

#include "framework/lua.h"

namespace opensnlua
{
/**Generic lua routine for the creation of solvers.
 * \param params ParameterBlock. A single block with at least one field
 *                   \"type\", which contains a registered solver type.
 * ## _
 *
 * Example:
 * \code
 * chiSolverCreate({type=cfem_diffusion.Solver})
 * \endcode
 *
 * \ingroup doc_PhysicsSolver
 * \deprecated This function is deprecated and will be removed soon.
 */
int chiSolverCreate(lua_State* L);

/** Initializes the solver at the given handle.
 *
 * \param solver_handle int Handle to the solver.
 *
 * \ingroup doc_PhysicsSolver
 * \author Jan
 */
int chiSolverInitialize(lua_State* L);

/** Executes the solver at the given handle.
 *
 * \param solver_handle int Handle to the solver.
 *
 * \ingroup doc_PhysicsSolver
 * \author Jan
 */
int chiSolverExecute(lua_State* L);

/** Performs a single timestep for the solver at the given handle.
 *
 * \param solver_handle int Handle to the solver.
 *
 * \ingroup doc_PhysicsSolver
 * \author Jan
 */
int chiSolverStep(lua_State* L);

/** Advances the time values of the solver at the given handle.
 *
 * \param solver_handle int Handle to the solver.
 *
 * \ingroup doc_PhysicsSolver
 * \author Jan
 */
int chiSolverAdvance(lua_State* L);

/** Sets a basic option of a solver.
 *
 * \param solver_handle int Handle to the reference solver.
 * \param option_name   string String-name of the option.
 * \param option_value  varying The value to assign to the option.
 *
 * \ingroup doc_PhysicsSolver
 * \deprecated This function is deprecated and will be removed soon.
 * \author Jan
 */
int chiSolverSetBasicOption(lua_State* L);

/** Returns the text name of the solver.
 *
 * \param solver_handle int Handle to the solver.
 *
 * \ingroup doc_PhysicsSolver
 * \author Jan
 */
int chiSolverGetName(lua_State* L);

/**Obtains a named list of the field functions associated with a solver.
 *
 * \param SolverHandle int A handle to the reference solver.
 *
 * \ingroup doc_PhysicsSolver
 */
int chiSolverGetFieldFunctionList(lua_State* L);

/** Returns arbitrary info specific for each solver.
 *
 * \param solver_handle int Handle to the solver.
 * \param info varying A single string or a table of values to call the solver
 * with.
 *
 * \ingroup doc_PhysicsSolver
 * \author Jan
 */
int chiSolverGetInfo(lua_State* L);

/**Sets a property of a solver.
 * \param handle int Solver handle.
 * \param property_table Table Table of properties to set. See solver specific
 * documentation.
 *
 * \ingroup doc_PhysicsSolver
 */
int chiSolverSetProperties(lua_State* L);
} // namespace opensnlua
