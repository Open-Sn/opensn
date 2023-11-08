#pragma once

#include "framework/chi_lua.h"

/** \defgroup LuaSurfaceMesh Surface Meshes
 * \ingroup LuaMesh
 *
 * Creates a new empty surface mesh.
 *
 * ### Example
 * Example usage:
 * \code
 * surfmesh = chiSurfaceMeshCreate()
 * \endcode
 *
 * \return Handle int Handle to the created surface mesh.
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int chiSurfaceMeshCreate(lua_State* L);

/** Loads mesh data from a wavefront object.
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
 * surfmesh1 = chiSurfaceMeshCreate()
 * chiSurfaceMeshImportFromOBJFile(surfmesh1, "MeshFile1.obj")
 *
 * -- Surface mesh used as Logical volume
 * lv_surfmesh1 = chiSurfaceMeshCreate()
 * chiSurfaceMeshImportFromOBJFile(lv_surfmesh1, "MeshFile3D.obj", false)
 *
 * lv1 = chiLogicalVolumeCreate(SURFACE, lv_surfmesh1)
 *
 * -- Surface mesh with transform
 * dx = 1.5
 * dy = -2.5
 * lv_surfmesh2 = chiSurfaceMeshCreate()
 * chiSurfaceMeshImportFromOBJFile(lv_surfmesh2, "MeshFile3D.obj", false,
 * {dx,dy,0.0})
 *
 * lv2 = chiLogicalVolumeCreate(SURFACE, lv_surfmesh2)
 * \endcode
 *
 * \return success bool Return true if file was successfully loaded and false
 *  otherwise.
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int chiSurfaceMeshImportFromOBJFile(lua_State* L);

/** Loads mesh data from a wavefront object.
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
int chiSurfaceMeshImportFromTriangleFiles(lua_State* L);

/** Exports mesh as a .obj format.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be
 * performed. \param FileName char* Path to the file to be exported.
 *
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int chiSurfaceMeshExportToObj(lua_State* L);

/** Exports mesh as a .poly format.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be
 * performed. \param FileName char* Path and basename to the file to be exported.
 *
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int chiSurfaceMeshExportPolyFile(lua_State* L);

/** Builds sweep ordering for a number of angles and checks whether any
 * cyclic dependencies are encountered.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be
 * performed. \param NumAngles int Number of azimuthal angles to use for checking
 * cycles.
 *
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int chiSurfaceMeshCheckCycles(lua_State* L);

/** Computes load balancing parameters for given predictive x and y cuts
 * without actually performing cuts.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be
 * performed. \param Xcuts table Array of x-values associated with the xcuts.
 * \param Ycuts table Array of y-values associated with the ycuts.
 *
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int chiComputeLoadBalancing(lua_State* L);
