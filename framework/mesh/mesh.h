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
class FieldFunctionInterpolationSlice;
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

enum MeshAttributes : int
{
  NONE = 0,
  ORTHOGONAL = (1 << 3),
  EXTRUDED = (1 << 4),
  UNSTRUCTURED = (1 << 5)
};

inline MeshAttributes
operator|(const MeshAttributes f1, const MeshAttributes f2)
{
  return static_cast<MeshAttributes>(static_cast<int>(f1) | static_cast<int>(f2));
}

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
