// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua
{

/**
 * Obtains a list of field functions, related only to scalar flux, from the transport solver.
 *
 * \param SolverIndex int Handle to the solver for which the list is to be obtained.
 *
 * \return Pair Table and count. Returns an array of handles and the amount of
 * elements in it (indexed from 1). \ingroup LBSLuaFunctions \author Jan
 */
int LBSGetScalarFieldFunctionList(lua_State* L);

/**
 * Obtain the power field function handle.
 *
 * \param SolverIndex int Handle to the solver for which the list is to be obtained.
 *
 * \return ff_handle The power field function handle.
 */
int LBSGetPowerFieldFunction(lua_State* L);

/**
 * Writes the angular fluxes of a LBS groupset to file.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 */
int LBSWriteGroupsetAngularFlux(lua_State* L);

/**
 * Reads the angular fluxes of a LBS groupset from a file.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param GroupsetIndex int Index to the groupset to which this function should
 *                          apply
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 */
int LBSReadGroupsetAngularFlux(lua_State* L);

/**
 * Writes the flux-moments of a LBS solution to file.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".h5"
 */
int LBSWriteFluxMoments(lua_State* L);

/**
 * Creates scattered source-moments, based on a LBS solution, and writes them
 * to file.
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".data"
 */
int LBSCreateAndWriteSourceMoments(lua_State* L);

/**
 * Reads flux-moments from a file and creates a scattering source from these
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
 */
int LBSReadFluxMomentsAndMakeSourceMoments(lua_State* L);

/**
 * Reads the source-moments from a file to a specific ext_src_moments_local-vector to be used
 * instead of a regular material/boundary source.
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
int LBSReadSourceMoments(lua_State* L);

/**
 * Reads flux-moments from a file to phi_old_local (the initial flux solution).
 *
 * \param SolverIndex int Handle to the solver for which the group
 * is to be created.
 *
 * \param file_base string Path+Filename_base to use for the output. Each location
 *                         will append its id to the back plus an extension ".h5"
 *
 * \param single_file_flag bool (Optional) Flag indicating that the file is a
 *                              single stand-alone file. The file_base will then
 *                              be used without adding the location-id, but still
 *                              with the ".data" appended. Default: false.
 */
int LBSReadFluxMoments(lua_State* L);

/**
 * Computes and returns the fission rate.
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
int LBSComputeFissionRate(lua_State* L);

/**
 * Initializes or reinitializes the materials. This normally happens
 * automatically during solver initialization but if the user wants to
 * swap/change XSs during the run then this will allow the material structures
 * to now deal with the new/changed materials.
 *
 * \param SolverIndex int Handle to the solver maintaining the information.
 *
 * \ingroup LBSLuaFunctions
 * \author Jan
 */
int LBSInitializeMaterials(lua_State* L);

} // namespace opensnlua
