#pragma once

/** Creates a 1D Mesh from an array of 1D vertices.
 *
 * \param x_nodes array_float An Array of floating point numbers denoting
 *                            1D nodes along x-axis.
 *
 * \return Two handles: unpartitioned-mesh, region
 *
 * \ingroup LuaMeshMacros
 *
 * ##_
 *
 * ### Example
 * An example 1D mesh creation below:
 * \code
 * chiMeshHandlerCreate()
 * nodes={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
 * chiMeshCreateUnpartitioned1DOrthoMesh(nodes)
 * chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
 * chiVolumeMesherExecute();
 * \endcode
 *
 * \author Nak
 */
int chiMeshCreateUnpartitioned1DOrthoMesh(lua_State* L);

/** Creates a 2D Orthogonal Mesh from arrays of 1D vertices.
 *
 * \param x_nodes array_float An Array of floating point numbers denoting
 *                            1D nodes along x-axis.
 * \param y_nodes array_float An Array of floating point numbers denoting
 *                            1D nodes along y-axis.
 *
 * \return Two handles: unpartitioned-mesh, region
 *
 * \ingroup LuaMeshMacros
 *
 * ##_
 *
 * ### Example
 * An example 2D mesh creation below:
 * \code
 * chiMeshHandlerCreate()
 * nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
 * nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
 * chiMeshCreateUnpartitioned2DOrthoMesh(nodesx,nodesy)
 * chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
 * chiVolumeMesherExecute();
 * \endcode
 *
 *  \author Nak
 */
int chiMeshCreateUnpartitioned2DOrthoMesh(lua_State* L);

/** Creates a 3D Orthogonal Mesh from arrays of 1D vertices. The
 * underlying mesher is an extruder.
 *
 * \param x_nodes array_float An Array of floating point numbers denoting
 *                            1D nodes along x-axis.
 * \param y_nodes array_float An Array of floating point numbers denoting
 *                            1D nodes along y-axis.
 * \param z_nodes array_float An Array of floating point numbers denoting
 *                            1D nodes along z-axis.
 *
 * \return Two handles: unpartitioned-mesh, region
 *
 * \ingroup LuaMeshMacros
 *
 * ##_
 *
 * ### Example
 * An example 3D mesh creation below:
 * \code
 * chiMeshHandlerCreate()
 * nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
 * nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
 * nodesz={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
 * chiMeshCreateUnpartitioned3DOrthoMesh(nodesx,nodesy,nodesz)
 * chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
 * chiVolumeMesherExecute();
 * \endcode
 *
 *  \author Nak
 */
int chiMeshCreateUnpartitioned3DOrthoMesh(lua_State* L);
