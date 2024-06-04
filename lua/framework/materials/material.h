// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua
{

/**
 * Adds a material to the problem. Materials are added to the global
 * physics handler and is therefore accessible across all meshes and solvers.
 *
 * \param Name char (Optional) Material name.
 *
 * \return MaterialHandle int Handle to the created material.
 *
 *
 * ##_
 *
 * ### Example\n
 * Example lua code:
 * \code
 * materials = {}
 * materials[1] = mat.AddMaterial("Test Material");
 * \endcode
 *
 *
 * \ingroup LuaPhysicsMaterials
 * \author Jan
 */
int MatAddMaterial(lua_State* L);

/**
 * Sets a material property for a given material.
 *
 * \param MaterialHandle int Index to the reference material.
 * \param PropertyType int Property type.
 * \param OperationType int Method used for setting the material property.
 * \param Information varying Varying information depending on the operation.
 *
 * ##_
 *
 * ###PropertyIndex\n
 * SCALAR_VALUE         =  Basic Scalar value property.\n
 * TRANSPORT_XSECTIONS  =  Multi-group transport cross section supporting numerous
 *                         features.\n
 * ISOTROPIC_MG_SOURCE = Isotropic Multigroup Source.\n
 *
 * ###OperationIndex\n
 * SINGLE_VALUE\n
 * Sets the property based on a single value. Requires a single value as additional
 * information. As a simple example consider the case where the user would like
 * to set a single constant thermal conductivity. This can be achieved with \n
 * FROM_ARRAY\n
 * Sets a property based on a Lua array indexed from 1 to N. Internally
 * will be converted to 0 to N-1. This method can be used to set mutligroup
 * cross sections or sources.
 * \n
 * SIMPLE_ONE_GROUP\n
 * Makes a simple, one-group material with a user-selectable amount of scattering.
 *  Expects two values:\n
 *  - float \f$\sigma_t \f$,
 *  - float scattering to total ratio (\f$c \f$)
 *
 * ####_
 *
 * OPENSN_XSFILE\n
 * Loads transport cross sections from OpenSn cross-section files. Expects
 * to be followed by a filepath specifying the xs-file.
 *
 * ####_
 *
 * OPENMC_XSLIB\n
 * Loads transport cross sections from an OpenMC MGXS cross-section library.
 * Expects to be followed by a string specifying the full path to the
 * cross-section library, a double representing the evaluation temperature
 * in K, and an optional string specifying the name of the cross-section
 * dataset.
 *
 * ####_
 *
 * EXISTING\n
 * Supply handle to an existing cross section and simply swap them out.
 *
 * \code
 * mat.SetProperty(materials[1],
 *                 TRANSPORT_XSECTIONS,
 *                 OPENSN_XSFILE,
 *                 "xs_3_170.xs",
 *                 "2518")
 * \endcode
 *
 * ##_
 *
 * ### Example 1
 * Simple scalar-valued property:
 * \code
 * materials = {}
 * materials[1] = mat.AddMaterial("Test Material");
 * mat.SetProperty(materials[0], SCALAR_VALUE, SINGLE_VALUE, 13.7)
 * \endcode
 *
 * where the thermal conductivity has been set to 13.7.\n
 *
 * ### Example 2
 * Isotropic multi-group source set from a lua table/array (12 groups):
 * \code
 * materials = {}
 * materials[1] = mat.AddMaterial("Test Material");
 *
 * num_groups = 12
 * src={}
 * for g=1,num_groups do
 *     src[g] = 0.0
 * end
 * mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)
 * \endcode
 *
 * ### Developer Info
 * Checklist for adding a new material property:
 *  - Under the "If user supplied name then find property index"-section
 *    add the appropriate code for setting the property index.
 *  - Add an else-if block for your property similar to the others. It should be
 *    intuitive if you look at the others.
 *  - Remember to add the filtering section if you need to support multiple type
 *    properties.
 *
 * \ingroup LuaPhysicsMaterials
 * \author Jan
 */
int MatSetProperty(lua_State* L);

/**
 * Returns a rich lua data-structure of the required property.
 *
 * \param MaterialHandle int Index to the reference material.
 * \param PropertyIndex int Property index. Or name of property.
 *
 *
 * \ingroup LuaPhysicsMaterials
 * \return Lua table of the desired property.
 */
int MatGetProperty(lua_State* L);

} // namespace opensnlua
