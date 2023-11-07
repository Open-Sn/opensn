#pragma once

#include "framework/chi_lua.h"

/**Creates a stand-alone transport cross section.
 *
 * \code
 * xs_graphite_clean = chiPhysicsTransportXSCreate()
 * chiPhysicsTransportXSSet(xs_grph_clean,
 *                          CHI_XSFILE,
 *                          "test/xs_graphite_pure.cxs")
 *
 * chiPhysicsMaterialSetProperty(materials[2],
 *                               TRANSPORT_XSECTIONS,
 *                               EXISTING,
 *                               xs_graphite_clean)
 * \endcode
 * \return Returns a handle to the cross section.
 *
 * \ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSCreate(lua_State* L);

/**Sets the properties of a transport cross section.
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
 * will be converted to 0 to N-1. This method can be used to set mutli-group
 * cross sections or sources.
 * \n
 * SIMPLEXS0\n
 * Makes a simple material with no transfer matrix just \f$\sigma_t \f$. Expects
 * two values: \n
 *  - int number of groups \f$G \f$,
 *  - float \f$\sigma_t \f$.
 *
 * ####_
 *
 * SIMPLEXS1\n
 * Makes a simple material with isotropic transfer matrix (L=0)
 * and mostly down scattering but with a few of the last groups
 * subject to up-scattering. Expects three values
 * values: \n
 *  - int number of groups (\f$G \f$),
 *  - float \f$\sigma_t \f$,
 *  - float scattering to total ratio (\f$c \f$)
 *
 * ####_
 *
 * CHI_XSFILE\n
 * Loads transport cross sections from CHI type cross section files. Expects
 * to be followed by a filepath specifying the xs-file.
 *
 *
 * ##_
 * ### Example
 * Example lua code:
 * \code
 * graphite = chiPhysicsTransportXSCreate()
 * chiPhysicsTransportXSSet(graphite,"xs_3_170.data","2518")
 * \endcode
 *
 *
 * \ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSSet(lua_State* L);

/**Makes a combined cross section from multiple other cross sections.
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
 * xs_1 = chiPhysicsTransportXSCreate()
 * xs_2 = chiPhysicsTransportXSCreate()
 * xs_3 = chiPhysicsTransportXSCreate()
 *
 * chiPhysicsTransportXSSet(xs_1,CHI_XSFILE,"test/xs_graphite_pure.cxs")
 * chiPhysicsTransportXSSet(xs_2,CHI_XSFILE,"test/xs_3_170.cxs")
 * chiPhysicsTransportXSSet(xs_3,CHI_XSFILE,"test/xs_air50RH.cxs")
 *
 * combo ={{xs_1, 0.5e5},
 *         {xs_2, 0.4e3},
 *         {xs_3, 0.3e2}}
 * aerated_graphite = chiPhysicsTransportXSMakeCombined(combo)
 *
 * chiPhysicsMaterialSetProperty(materials[1],
 *                               TRANSPORT_XSECTIONS,
 *                               EXISTING,
 *                               aerated_graphite)
 * \endcode
 *
 *  \return Returns a handle to another cross section object that contains the
 *          desired combination.
 *
 * \ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSMakeCombined(lua_State* L);

/**Sets a combined cross section from multiple other cross sections. This
 *  function can be called multiple times on the same cross section handle.
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
 * xs_1 = chiPhysicsTransportXSCreate()
 * xs_2 = chiPhysicsTransportXSCreate()
 * xs_3 = chiPhysicsTransportXSCreate()
 *
 * chiPhysicsTransportXSSet(xs_1,CHI_XSFILE,"test/xs_graphite_pure.cxs")
 * chiPhysicsTransportXSSet(xs_2,CHI_XSFILE,"test/xs_3_170.cxs")
 * chiPhysicsTransportXSSet(xs_3,CHI_XSFILE,"test/xs_air50RH.cxs")
 *
 * combo ={{xs_1, 0.5e5},
 *         {xs_2, 0.4e3},
 *         {xs_3, 0.3e2}}
 * aerated_graphite = chiPhysicsTransportXSMakeCombined(combo)
 *
 * chiPhysicsTransportXSSetCombined(aerated_graphite,combo)
 * \endcode
 *
 * \ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSSetCombined(lua_State* L);

/**Obtains a lua table of all the cross section values.
 *
 * \param XS_handle int Handle to the cross section to be modified.
 *
 * ## _
 *
 * To print the contents of the table, execute the following:
 * \code
 * xs = chiPhysicsTransportXSGet(xs_handle)
 * for i,v in pairs(xs) do
 *     print(i,v)
 * end
 * \endcode
 * \ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSGet(lua_State* L);

/** Exports a cross section to ChiTech format.
 *
 * \param XS_handle int Handle to the cross section to be exported.
 * \param file_name string The name of the file to which the XS is to be exported.
 *
 * \ingroup LuaTransportXSs
 */
int chiPhysicsTransportXSExportToChiTechFormat(lua_State* L);
