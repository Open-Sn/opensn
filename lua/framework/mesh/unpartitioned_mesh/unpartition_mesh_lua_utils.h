// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "lua/framework/lua.h"

namespace opensnlua
{

/**
 * Creates an unpartitioned mesh from VTK Unstructured mesh files.
 *
 * \param file_name char Filename of the .vtu file.
 * \param field char Name of the cell data field from which to read
 *                   material and boundary identifiers (optional).
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = mesh.UnpartitionedMeshFromVTU("ZMeshTest_0.vtu")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
 * \endcode
 *
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int MeshUnpartitionedMeshFromVTU(lua_State* L);

/**
 * Creates an unpartitioned mesh from VTK Partitioned Unstructured mesh files
 * (.pvtu).
 *
 * \param file_name char Filename of the .vtu file.
 * \param field char Name of the cell data field from which to read
 *                   material and boundary identifiers (optional).
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = mesh.UnpartitionedMeshFromPVTU("ZMeshTest_0.vtu")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
 * \endcode
 *
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int MeshUnpartitionedMeshFromPVTU(lua_State* L);

/**
 * Creates an unpartitioned mesh from starccm+ exported Ensight Gold mesh files.
 *
 * \param file_name char Filename of the .case file.
 * \param scale float Scale to apply to the mesh
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = mesh.UnpartitionedMeshFromEnsightGold("resources/TestObjects/Sphere.case")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
 * \endcode
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int MeshUnpartitionedMeshFromEnsightGold(lua_State* L);

/**
 * Creates an unpartitioned mesh from a wavefront .obj file.
 *
 * \param file_name char Filename of the .case file.
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = mesh.UnpartitionedMeshFromWavefrontOBJ("resources/TestObjects/TriangleMesh2x2.obj")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
 * \endcode
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int MeshUnpartitionedMeshFromWavefrontOBJ(lua_State* L);

/**
 * Creates an unpartitioned mesh from a .msh file.
 *
 * \param file_name char Filename of the .msh file.
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = mesh.UnpartitionedMeshFromMshFormat("File.msh")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
 * \endcode
 *
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int MeshUnpartitionedMeshFromMshFormat(lua_State* L);

/**
 * Creates an unpartitioned mesh from ExodusII format.
 *
 * \param file_name char Filename of the .case file.
 * \param scale float Scale to apply to the mesh
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = mesh.UnpartitionedMeshFromExodusII("resources/TestObjects/Mesh.e")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
 * \endcode
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int MeshUnpartitionedMeshFromExodusII(lua_State* L);

} // namespace opensnlua
