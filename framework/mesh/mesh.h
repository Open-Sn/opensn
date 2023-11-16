#pragma once

#include <vector>
#include <iostream>
#include <memory>

/**
 * Meshes in ChiTech follow the concept of Regions. In any given region the
 * boundaries are a collection of either line-meshes (2D) or
 * surface-meshes (3D).
 */
namespace opensn
{

struct Vector3;
typedef Vector3 Normal;
typedef Vector3 Vertex;

struct Matrix3x3;

struct Face;
struct Edge;
struct PolyFace;

class SPDS;

// Cells
class Cell;

// Field function interpolation
class FieldFunctionInterpolation;
class FieldFunctionInterpolationSlice;
class FieldFunctionInterpolationLine;
class FieldFunctionInterpolationVolume;

// Meshes
class SurfaceMesh;
class UnpartitionedMesh;
class MeshContinuum;
typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
typedef std::shared_ptr<const MeshContinuum> MeshContinuumConstPtr;

// Logical Volumes
class LogicalVolume;
class SphereLogicalVolume;
class RPPLogicalVolume;
class RCCLogicalVolume;
class SurfaceMeshLogicalVolume;
class BooleanLogicalVolume;

// Mesh handler
class MeshHandler;

// Surface Meshers
class SurfaceMesher;
class SurfaceMesherPassthrough;
class SurfaceMesherPredefined;

// Volume meshers
class VolumeMesher;
class VolumeMesherExtruder;
class VolumeMesherPredefinedUnpartitioned;

enum MeshAttributes : int
{
  NONE = 0,
  DIMENSION_1 = (1 << 0),
  DIMENSION_2 = (1 << 1),
  DIMENSION_3 = (1 << 2),
  ORTHOGONAL = (1 << 3),
  EXTRUDED = (1 << 4),
  UNSTRUCTURED = (1 << 5)
};

inline MeshAttributes
operator|(const MeshAttributes f1, const MeshAttributes f2)
{
  return static_cast<MeshAttributes>(static_cast<int>(f1) | static_cast<int>(f2));
}

/**
 * Obtains a reference to the current mesh handler from the global stack.
 *
 * If the stack is empty this routine will through `std::logic_error`.
 * \author Jan
 */
MeshHandler& GetCurrentHandler();

/**
 * Adds a new mesh handler to the stack, sets it as the current handler and returns a handle to it.
 */
size_t PushNewHandlerAndGetIndex();

// Domain decompositions

/**Makes a centroid based load balance factor calculation.
 *
 * \author Jan*/
double
ComputeLBF(std::vector<Vector3>& points, std::vector<double>& x_cuts, std::vector<double>& y_cuts);

/** Decomposes a 2D surface mesh using the centroids in a Px-Py fashion.*/
void DecomposeSurfaceMeshPxPy(const SurfaceMesh& smesh, int Px, int Py);

/**
 * Creates a 1D slab mesh from a set of vertices.
 */
size_t CreateUnpartitioned1DOrthoMesh(std::vector<double>& vertices_1d);

/**
 * Creates a 2D orthogonal mesh from a set of vertices in x and y.
 * The vertices along a dimension merely represents the divisions. They are not the complete
 * vertices defining a cell. For example:
 * \code
 * std::vector<Vertex> vertices_x = {0.0,1.0,2.0};
 * std::vector<Vertex> vertices_y = {0.0,1.0,2.0};
 * opensn::CreateUnpartitioned2DOrthoMesh(vertices_x,vertices_y);
 * \endcode
 *
 * This code will create a 2x2 mesh with \f$ \vec{x} \in [0,2]^2 \f$.
 *
 */
size_t CreateUnpartitioned2DOrthoMesh(std::vector<double>& vertices_1d_x,
                                      std::vector<double>& vertices_1d_y);

/**
 * Creates a 3D orthogonal mesh from a set of vertices in x,y,z.
 * The vertices along a dimension merely represents the divisions. They are not the complete
 * vertices defining a cell. For example:
 * \code
 * std::vector<Vertex> vertices_x = {0.0,1.0,2.0};
 * std::vector<Vertex> vertices_y = {0.0,1.0,2.0};
 * std::vector<Vertex> vertices_z = {0.0,1.0,2.0};
 * opensn::CreateUnpartitioned3DOrthoMesh(vertices_x,vertices_y,vertices_z);
 * \endcode
 *
 * This code will create a 2x2 mesh with \f$ \vec{x} \in [0,2]^2 \f$.
 *
 */
size_t CreateUnpartitioned3DOrthoMesh(std::vector<double>& vertices_1d_x,
                                      std::vector<double>& vertices_1d_y,
                                      std::vector<double>& vertices_1d_z);
} // namespace opensn

#include "framework/mesh/mesh_vector.h"
#include "framework/mesh/mesh_matrix3x3.h"
#include "framework/mesh/mesh_face.h"
#include "framework/mesh/mesh_edge_loops.h"

#include "framework/mesh/sweep_utilities/sweep_namespace.h"
