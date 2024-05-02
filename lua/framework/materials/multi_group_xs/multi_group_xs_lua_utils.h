// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/lua.h"

namespace opensnlua
{

/**
 * Creates a stand-alone transport cross section.
 *
 * \code
 * xs_graphite_clean = xs.Create()
 * xs.Set(xs_grph_clean, OPENSN_XSFILE, "test/xs_graphite_pure.xs")
 *
 * mat.SetProperty(materials[2],
 *                               TRANSPORT_XSECTIONS,
 *                               EXISTING,
 *                               xs_graphite_clean)
 * \endcode
 * \return Returns a handle to the cross section.
 *
 * \ingroup LuaTransportXSs
 */
int XSCreate(lua_State* L);

/**
 * Sets the properties of a transport cross section.
 *
 *  \param XS_handle int Handle to the cross section to be modified.
 *  \param OperationIndex int Method used for setting the xs property.
 *  \param Information varying Varying information depending on the operation.
 *
 *  ##_
 *
 * ###OperationIndex
 * SINGLE_VALUE\n
 * Sets the property based on a single value. Requires a single value as additional
 * information. As a simple example consider the case where the user would like
 * to set a single constant thermal conductivity. This can be achieved with \n
 * FROM_ARRAY\n
 * Sets a property based on a Lua array indexed from 1 to N. Internally
 * will be converted to 0 to N-1. This method can be used to set multi-group
 * cross sections or sources.
 * \n
 * SIMPLE_ONE_GROUP\n
 * Makes a simple, one-group material with a user-selectable amount of scattering.
 *  Expects two values: \n
 *  - float \f$\sigma_t \f$,
 *  - float scattering to total ratio (\f$c \f$)
 *
 * ####_
 *
 * OPENSN_XSFILE\n
 * Loads transport cross sections from OpenSn cross-section files. Expects
 * to be followed by a filepath specifying the xs-file.
 *
 *
 * ##_
 * ### Example
 * Example lua code:
 * \code
 * graphite = xs.Create()
 * xs.Set(graphite, "xs_3_170.data", "2518")
 * \endcode
 *
 *
 * \ingroup LuaTransportXSs
 */
int XSSet(lua_State* L);

/**
 * Makes a combined cross section from multiple other cross sections.
 *
 *  \param Combinations table A lua-table with each element another table
 *                            containing a handle to an existing xs and a
 *                            scalar multiplier.
 *
 *  ## _
 *
 * ###Example:
 * Example lua code:
 * \code
 * xs_1 = xs.Create()
 * xs_2 = xs.Create()
 * xs_3 = xs.Create()
 *
 * xs.Set(xs_1, OPENSN_XSFILE, "test/xs_graphite_pure.xs")
 * xs.Set(xs_2, OPENSN_XSFILE, "test/xs_3_170.xs")
 * xs.Set(xs_3, OPENSN_XSFILE, "test/xs_air50RH.xs")
 *
 * combo ={{xs_1, 0.5e5},
 *         {xs_2, 0.4e3},
 *         {xs_3, 0.3e2}}
 * aerated_graphite = xs.MakeCombined(combo)
 *
 * mat.SetProperty(materials[1],
 *                               TRANSPORT_XSECTIONS,
 *                               EXISTING,
 *                               aerated_graphite)
 * \endcode
 *
 *  \return Returns a handle to another cross-section object that contains the
 *          desired combination.
 *
 * \ingroup LuaTransportXSs
 */
int XSMakeCombined(lua_State* L);

/**
 * Sets a combined cross section from multiple other cross sections. This function can be called
 * multiple times on the same cross-section handle.
 *
 *  \param XS_handle int Handle to the cross section to be modified.
 *  \param Combinations table A lua-table with each element another table
 *                            containing a handle to an existing xs and a
 *                            scalar multiplier.
 *
 *  ## _
 *
 * ###Example:
 * Example lua code:
 * \code
 * xs_1 = xs.Create()
 * xs_2 = xs.Create()
 * xs_3 = xs.Create()
 *
 * xs.Set(xs_1, OPENSN_XSFILE, "test/xs_graphite_pure.xs")
 * xs.Set(xs_2, OPENSN_XSFILE, "test/xs_3_170.xs")
 * xs.Set(xs_3, OPENSN_XSFILE, "test/xs_air50RH.xs")
 *
 * combo ={{xs_1, 0.5e5},
 *         {xs_2, 0.4e3},
 *         {xs_3, 0.3e2}}
 * aerated_graphite = xs.MakeCombined(combo)
 *
 * xs.SetCombined(aerated_graphite,combo)
 * \endcode
 *
 * \ingroup LuaTransportXSs
 */
int XSSetCombined(lua_State* L);

/**
 * Obtains a lua table of all the cross-section values.
 *
 * \param XS_handle int Handle to the cross section to be modified.
 *
 * ## _
 *
 * To print the contents of the table, execute the following:
 * \code
 * xs = xs.Get(xs_handle)
 * for i,v in pairs(xs) do
 *     print(i,v)
 * end
 * \endcode
 * \ingroup LuaTransportXSs
 */
int XSGet(lua_State* L);

/**
 * Exports a cross section to OpenSn format.
 *
 * \param XS_handle int Handle to the cross section to be exported.
 * \param file_name string The name of the file to which the XS is to be exported.
 *
 * \ingroup LuaTransportXSs
 */
int XSExportToOpenSnFormat(lua_State* L);

} // namespace opensnlua
