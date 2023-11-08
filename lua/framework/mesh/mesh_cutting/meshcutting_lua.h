#pragma once

#include "framework/chi_lua.h"

/**Cuts a mesh.
 *
 * \param PlanePoint Vec3 A 3 element lua-table containing the reference
 *                        point of the plane.
 * \param PlaneNormal Vec3 A 3 element lua-table containing the plane
 *                         normal. Will be normalized when applied.
 * \param MergeTol double (Optional). Vertices within a distance MergeTol from
 *                        the plane will be snapped to the plane. Default: 1.0e-3
 * \param GeneralTol double (Optional). A general tolerance to be used for
 *                          comparing whether two real values are equal.
 *                          Default: 1.0e-10.
 */
int chiCutMesh(lua_State* L);

/**Counts the number of cells with a logical volume.*/
int chiCountMeshInLogicalVolume(lua_State* L);
