#pragma once

#include <vector>
#include <iostream>
#include <memory>

/**
 * Namespace for all meshing features
 *
 * Meshes in ChiTech follow the concept of Regions. In any given region the
 * boundaries are a collection of either line-meshes (2D) or
 * surface-meshes (3D).
 */
namespace chi_mesh
{
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FORWARD DECLARATIONS
struct Vector3;
typedef Vector3 Normal;
typedef Vector3 Vertex;

struct Matrix3x3;
struct TensorRank2Dim3;

struct Face;
struct Edge;
struct PolyFace;

struct SPDS;

//=================================== Cells
class Cell;

//=================================== Field function interpolation
class FieldFunctionInterpolation;
class FieldFunctionInterpolationSlice;
class FieldFunctionInterpolationLine;
class FieldFunctionInterpolationVolume;

//=================================== Meshes
class SurfaceMesh;
class UnpartitionedMesh;
class MeshContinuum;
typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
typedef std::shared_ptr<const MeshContinuum> MeshContinuumConstPtr;

//=================================== Logical Volumes
class LogicalVolume;
class SphereLogicalVolume;
class RPPLogicalVolume;
class RCCLogicalVolume;
class SurfaceMeshLogicalVolume;
class BooleanLogicalVolume;

//=================================== Mesh handler
class MeshHandler;

//=================================== Surface Meshers
class SurfaceMesher;
class SurfaceMesherPassthrough;
class SurfaceMesherPredefined;

//==================================== Volume meshers
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

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ROUTINES

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

//=================================== Domain decompositions
double
ComputeLBF(std::vector<Vector3>& points, std::vector<double>& x_cuts, std::vector<double>& y_cuts);
void DecomposeSurfaceMeshPxPy(const SurfaceMesh& smesh, int Px, int Py);

size_t CreateUnpartitioned1DOrthoMesh(std::vector<double>& vertices_1d);

size_t CreateUnpartitioned2DOrthoMesh(std::vector<double>& vertices_1d_x,
                                      std::vector<double>& vertices_1d_y);

size_t CreateUnpartitioned3DOrthoMesh(std::vector<double>& vertices_1d_x,
                                      std::vector<double>& vertices_1d_y,
                                      std::vector<double>& vertices_1d_z);
} // namespace chi_mesh

#include "chi_meshvector.h"
#include "chi_meshmatrix3x3.h"
#include "chi_meshtensor_rank2_dim3.h"
#include "chi_meshface.h"
#include "chi_mesh_edgeloops.h"

#include "mesh/SweepUtilities/sweep_namespace.h"
