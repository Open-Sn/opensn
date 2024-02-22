#pragma once

#include "framework/lua.h"

/**
 * Creates a logical volume.
 *
 * \param TypeIndex int Volume type.
 * \param Values varying Parameters.
 *
 * \warning This function is deprecated
 *
 * ##_
 *
 * ### TypeIndex
 * SPHERE_ORIGIN =  Sphere at the origin. [Requires: R]\n
 * SPHERE  =  Sphere at user supplied location. [Requires: x,y,z,R]\n
 * RPP=  Rectangular ParalleliPiped. [Requires: xmin,xmax,ymin,ymax,zmin,zmax]\n
 * RCC=  Right Circular Cylinder. [Requires: x0,y0,z0, vx,vy,vz and R]\n
 * SURFACE= Logical volume determined from a surface mesh. [Requires: a handle
 *     to a surface mesh]\n
 * BOOLEAN= Boolean combination of other volumes.
 *  [Requires pairs of values: bool,volumeHandle]\n
 *
 * ### Examples
 * Example usage:
 * \code
 * -- Sphere at origin radius 1.0
 * lv1 = logvol.Create(SPHERE_ORIGIN, 1.0)
 *
 * -- Sphere centered at (0.1,0.2,0.3) with radius 1.0
 * lv2 = logvol.Create(SPHERE, 0.1, 0.2, 0.3, 1.0)
 *
 * -- Rectangular parallelepiped
 * xmin = -1.0; xmax = 1.0
 * ymin = -2.0; ymax = 1.0
 * zmin = -1000.0, zmax = 1000.0
 * lv3 = logvol.Create(RPP, xmin, xmax, ymin, ymax, zmin, zmax)
 *
 * -- Right Circular Cylinder
 * basex = 1.0; basey = 0.0; basez = 0.5
 * upvecx = 0.0; upvecy = 0.0; upvecz = 1.0
 * R = 1.0
 * lv4 = logvol.Create(RPP, basex, basey, basez, upvecx, upvecy, upvecz,
 * R)
 *
 * -- Surface mesh
 * lv_surfmesh = SurfaceMeshCreate()
 * SurfaceMeshImportFromOBJFile(lv_surfmesh, "MeshFile3D.obj", false)
 *
 * lv5 = logvol.Create(SURFACE, lv_surfmesh)
 *
 * -- Boolean combination
 * lv6 = logvol.Create(BOOLEAN, {{true , lv5},  -- inside logical volume 5
 *                                        {false, lv1},  -- outside logical volume
 * 1 {false, lv2}}) -- outside logical volume 2 \endcode
 *
 * \return Handle int Handle to the created logical volume.
 * \ingroup LuaLogicVolumes
 * \author Jan
 */
int LogVolCreate(lua_State* L);

/**
 * Evaluates whether a point is within the logical volume.
 *
 * \param LVHandle int Handle to the logical volume.
 * \param Point_x double X-coordinate of the point.
 * \param Point_y double Y-coordinate of the point.
 * \param Point_z double Z-coordinate of the point.
 *
 * \return Sense true if inside the logical volume and false if outside.
 * \ingroup LuaLogicVolumes
 */
int LogVolPointSense(lua_State* L);
