#pragma once

#include "framework/lua.h"

/** Creates a surface preprocessor.
 *
 * \param SurfaceMesherType int Surface Remesher type. See SurfaceMesherType.
 *
 * ## _
 *
 * ###SurfaceMesherType:\n
 * SURFACEMESHER_PREDEFINED\n
 *  Makes no modification to the region surfaces.\n\n
 *
 * \code
 * chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * \endcode
 *
 * ## _
 *
 * \ingroup LuaSurfaceMesher
 * \author Jan
 */
int chiSurfaceMesherCreate(lua_State* L);

/** Executes the surface meshing pipeline.
 *
 * \ingroup LuaSurfaceMesher
 * \author Jan
 */
int chiSurfaceMesherExecute(lua_State* L);

/** Sets a property of a surface mesher.
 *
 * \param PropertyNumber int Handle of the property to be set.
 * \param PropertyValue varying Value of the property.
 *
 * Properties:\n
 *  MAX_AREA = Area constraint.\n
 *
 * \ingroup LuaSurfaceMesher
 * \author Jan
 */
int chiSurfaceMesherSetProperty(lua_State* L);
// int chiSurfaceMesherExportToObj(lua_State* L);
