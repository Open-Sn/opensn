// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua
{

/**
 * Creates a new empty surface mesh.
 *
 * ### Example
 * Example usage:
 * \code
 * surfmesh = mesh.SurfaceMeshCreate()
 * \endcode
 *
 * \return Handle int Handle to the created surface mesh.
 * \ingroup LuaSurfaceMesh
 * \defgroup LuaSurfaceMesh Surface Meshes
 * \author Jan
 */
int MeshSurfaceMeshCreate(lua_State* L);

/**
 * Loads mesh data from a wavefront object.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be
 * performed. \param FileName string Path to the file to be imported. \param
 * polyflag bool (Optional)Flag indicating whether triangles are to be read as
 * polygons. [Default: true (read as polygons)]. \param transform table3 (Optional)
 * Translation vector to move all the vertices. [Default: none].
 *
 * ### Note:
 * If the intent of a surface mesh is to serve as a 3D logical volume then
 * the `polyFlag` parameter should be set to false.
 *
 * ### Example
 * Example usage:
 * \code
 * -- Basic example
 * surfmesh1 = mesh.SurfaceMeshCreate()
 * mesh.SurfaceMeshImportFromOBJFile(surfmesh1, "MeshFile1.obj")
 *
 * -- Surface mesh used as Logical volume
 * lv_surfmesh1 = mesh.SurfaceMeshCreate()
 * mesh.SurfaceMeshImportFromOBJFile(lv_surfmesh1, "MeshFile3D.obj", false)
 *
 * lv1 = logvol.SurfaceMeshLogicalVolume.Create({surface_mesh_handle=lv_surfmesh1})
 *
 * -- Surface mesh with transform
 * dx = 1.5
 * dy = -2.5
 * lv_surfmesh2 = mesh.SurfaceMeshCreate()
 * mesh.SurfaceMeshImportFromOBJFile(lv_surfmesh2, "MeshFile3D.obj", false,
 * {dx,dy,0.0})
 *
 * lv2 = logvol.SurfaceMeshLogicalVolume.Create({surface_mesh_handle=lv_surfmesh2})
 *
 * \endcode
 *
 * \return success bool Return true if file was successfully loaded and false
 *  otherwise.
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int MeshSurfaceMeshImportFromOBJFile(lua_State* L);

/**
 * Loads mesh data from a wavefront object.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be
 * performed. \param FileName char* Path to the file to be imported. \param
 * polyflag bool (Optional)Flag indicating whether triangles are to be read as
 * polygons. [Default: true)
 *
 * \return success bool Return true if file was successfully loaded and false
 *  otherwise.
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int MeshSurfaceMeshImportFromTriangleFiles(lua_State* L);

/**
 * Computes load balancing parameters for given predictive x and y cuts without actually performing
 * cuts.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be performed.
 * \param Xcuts table Array of x-values associated with the xcuts.
 * \param Ycuts table Array of y-values associated with the ycuts.
 *
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int MeshComputeLoadBalancing(lua_State* L);

} // namespace opensnlua
