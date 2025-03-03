// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/vector3.h"
#include <array>

namespace opensn
{

/// Data structure for a triangular face.
struct Face
{
  /// vertex indices
  std::array<int, 3> v_index{-1, -1, -1};
  /// normal indices
  std::array<int, 3> n_index{-1, -1, -1};
  /// vertex texture indices
  std::array<int, 3> vt_index{-1, -1, -1};
  /// edge indices
  std::array<std::array<int, 4>, 3> e_index{{{-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}}};

  Vector3 geometric_normal;
  Vector3 assigned_normal;
  Vector3 face_centroid;

  bool invalidated{false};
};

/**
 * Data structure for a polygon face.
 *
 * edges
 * An array of 4 integers.
 * [0] = Vertex index of edge start.
 * [1] = Vertex index of edge end.
 * [2] = Index of the face adjoining this edge (not the current face).
 *       -1 if not connected to anything,-1*boundary_index if connected
 *       to a boundary.
 * [3] = Edge number of adjoining face. -1 if not connected
 *       to anything. 0 if a boundary.
 *
 * face_indices
 *  [0] = Index of the adjoining cell. -1 if not connected to anything.
 *        -1*boundary_index if connected to a boundary.
 *  [1] = Face number of adjoining cell. -1 if not connected
 *       to anything. 0 if a boundary.
 *  [2] = Partition ID of adjacent cell.
 */
struct PolyFace
{
  std::vector<int> v_indices;
  std::vector<std::vector<int>> edges;
  std::array<int, 3> face_indices{{-1, -1, -1}};

  Vector3 geometric_normal;
  Vector3 face_centroid;

  bool invalidated{false};
};

} // namespace opensn
