// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <iostream>
#include <memory>

/**
 * Meshes in OpenSn follow the concept of Regions. In any given region the
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
class FieldFunctionInterpolationLine;
class FieldFunctionInterpolationVolume;

// Meshes
class SurfaceMesh;
class UnpartitionedMesh;
class MeshContinuum;
// Logical Volumes
class LogicalVolume;
class SphereLogicalVolume;
class RPPLogicalVolume;
class RCCLogicalVolume;
class SurfaceMeshLogicalVolume;
class BooleanLogicalVolume;

// Volume meshers
class VolumeMesher;
class VolumeMesherExtruder;
class VolumeMesherPredefinedUnpartitioned;

enum MeshType : int
{
  ORTHOGONAL,
  UNSTRUCTURED
};

enum BoundaryID : int
{
  XMIN = 0,
  XMAX = 1,
  YMIN = 2,
  YMAX = 3,
  ZMIN = 4,
  ZMAX = 5
};

struct OrthoMeshAttributes
{
  size_t Nx = 0;
  size_t Ny = 0;
  size_t Nz = 0;
};

/**
 * Obtains the current mesh from the global stack.
 */
std::shared_ptr<MeshContinuum> GetCurrentMesh();

} // namespace opensn

#include "framework/mesh/mesh_vector.h"
#include "framework/mesh/mesh_matrix3x3.h"
#include "framework/mesh/mesh_face.h"
#include "framework/mesh/mesh_edge_loops.h"
