#pragma once

#include "framework/lua.h"

/** \defgroup LuaDomainDecomposition Domain decomposition
 * \ingroup LuaMesh
 *
 * Decomposes a surface mesh into block px py elements.
 * \image html "InProgressImage.png" width=200px
 *
 * \param Surface mesh handler
 * \param Px int Number of divisions in x.
 * \param Py int Number of divisions in y.
 *
 * \ingroup LuaDomainDecomposition
 * \author Jan
 */
int chiDecomposeSurfaceMeshPxPy(lua_State* L);
