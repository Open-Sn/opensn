// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/lua.h"

namespace opensnlua
{

/** Creates a new field function interpolation.
 *
 * \param FFITypeIndex int Type of field function interpolation.
 *
 * ##_
 *
 * ###FFITypeIndex\n
 * POINT           = A point probe. \n
 * SLICE           = Two dimensional slice of the mesh. \n
 * LINE            = Line defined by two points.\n
 * VOLUME          = Volume either referring to the entire volume or that of a
 *                   logical volume assigned to the interpolator.\n
 *
 * \return Handle int Handle to the created interpolation.
 * \ingroup LuaFFInterpol
 * \author Jan
 */
int FFInterpolationCreate(lua_State* L);

/** Creates a new field function interpolation.
 *
 * \param FFIHandle int Handle to the field function interpolation.
 * \param PropertyIndex int Type of property to set.
 *
 * ##_
 *
 * ###PropertyIndex\n
 * ADD_FIELDFUNCTION     = Add field function to interpolation.\n
 * SET_FIELDFUNCTIONS    = Sets the field functions to interpolate using a lua
 * table.\n PROBEPOINT            = Reference point for POINT type FFIs.\n
 * SLICE_POINT           = Reference point for SLICE type FFIs.\n
 * SLICE_NORMAL          = Normal of slice plane.\n
 * SLICE_TANGENT         = Tangent vector of slice plane.\n
 * SLICE_BINORM          = Binorm vector of slice plane.\n
 * LINE_FIRSTPOINT   = Line start point.\n
 * LINE_SECONDPOINT  = Line end point.\n
 * LINE_NUMBEROFPOINTS = Number of points to put in the line interpolator.
 *                           Minimum 2.\n
 * LINE_CUSTOM_ARRAY = Adds custom array to line interpolator.\n
 * OPERATION  =  Some interpolations support operation types. See OpTypes.\n
 * LOGICAL_VOLUME = To be followed by a handle to a logical volume to be
 *                  used by the interpolator.\n
 *
 * ###OpTypes
 * Basic operations are volume sum, volume average or volume max. The volume
 * sum is computed from
 * \f[
 * Sum = \sum_k \sum_i u_i \int_V N_i .dV.
 * \f]
 * The volume average is computed from
 * \f[
 * Avg = \frac{\sum_k \sum_i u_i \int_V N_i .dV}{\sum_k \sum_i \int_V N_i .dV}
 * \f]
 * The volume max is simply \f$ max_{k,i}(u_i) \f$.\n\n
 *
 * OP_SUM\n
 * For volume interpolations, computes the volume integral.\n
 * \n
 * OP_AVG\n
 * For volume interpolations, computes the volume average.\n
 * OP_MAX\n
 * For volume interpolations, computes the volume max.\n
 * \n
 * A modified version of these operations are also available. Instead of OP_SUM,
 * OP_AVG and OP_MAX, the user may supply OP_SUM_LUA, OP_AVG_LUA and OP_MAX_LUA
 * which then needs to be followed by a string value `LuaFunctionName` of a lua
 * function of the following form:
 *
 * \code
 * function LuaFunctionName(ff_value, mat_id)
 *  ret_val = 0.0;   --Or some computation
 *  return ret_val
 * end
 * \endcode
 *
 * This code will be called to return a value \f$ f(u_i) \f$ to be used instead of
 * the field function \f$ u_i \f$.\n
 *
 * Example:
 * \code
 * xwing=2.0
 * function IntegrateMaterialVolume(ff_value,mat_id)
 *     return xwing
 * end
 * ffi2 = fieldfunc.FFInterpolationCreate(VOLUME)
 * curffi = ffi2
 * fieldfunc.SetProperty(curffi,OPERATION,OP_SUM_LUA,"IntegrateMaterialVolume")
 * fieldfunc.SetProperty(curffi,LOGICAL_VOLUME,vol0)
 * fieldfunc.SetProperty(curffi,ADD_FIELDFUNCTION,fftemp)
 *
 * fieldfunc.Initialize(curffi)
 * fieldfunc.Execute(curffi)
 * print(fieldfunc.GetValue(curffi))
 * \endcode
 *
 * The code above will return 2.0 times the volume of cells included in the logical
 * volume `vol0`.
 *
 *
 * \return Handle int Handle to the created interpolation.
 * \ingroup LuaFFInterpol
 * \author Jan
 */
int FFInterpolationSetProperty(lua_State* L);

/** Initialize interpolator.
 *
 * \param FFIHandle int Handle to the field function interpolation.
 *
 * \ingroup LuaFFInterpol
 * \author Jan
 */
int FFInterpolationInitialize(lua_State* L);

/** Execute interpolator.
 *
 * \param FFIHandle int Handle to the field function interpolation.
 *
 * \ingroup LuaFFInterpol
 * \author Jan
 */
int FFInterpolationExecute(lua_State* L);

/** Export interpolation to python line,contour plot depending on the
 * type of interpolation.
 *
 * \param FFIHandle int Handle to the field function interpolation.
 * \param BaseName char Base name to be used for exported files.
 *
 * \ingroup LuaFFInterpol
 * \author Jan
 */
int FFInterpolationExportPython(lua_State* L);

/** Gets the value(s) associated with an interpolation provided the
 * interpolation type has an associated value.
 *
 * \param FFIHandle int Handle to the field function interpolation.
 *
 * ###Note:
 * Currently only the POINT, LINE and VOLUME interpolation supports obtaining a
 * value. For the POINT and VOLUME types a single value is returned. For the LINE
 * type a table of tables is returned with the first index being the field function
 * (in the order it was assigned) and the second index being the point index.
 *
 * \ingroup LuaFFInterpol
 * \author Jan
 */
int FFInterpolationGetValue(lua_State* L);

} // namespace opensnlua
