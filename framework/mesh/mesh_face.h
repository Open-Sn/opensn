// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

/// Data structure for a triangular face.
struct Face
{
  int v_index[3];
  int n_index[3];
  int vt_index[3];
  int e_index[3][4];

  Vector3 geometric_normal;
  Vector3 assigned_normal;
  Vector3 face_centroid;

  bool invalidated;

  Face()
  {
    for (int k = 0; k < 3; ++k)
    {
      v_index[k] = -1;
      n_index[k] = -1;
      vt_index[k] = -1;
      e_index[k][0] = -1;
      e_index[k][1] = -1;
      e_index[k][2] = -1;
      e_index[k][3] = -1;
      invalidated = false;
    }
  }

  void SetIndices(int a, int b, int c)
  {
    v_index[0] = a;
    v_index[1] = b;
    v_index[2] = c;

    e_index[0][0] = a;
    e_index[0][1] = b;
    e_index[1][0] = b;
    e_index[1][1] = c;
    e_index[2][0] = c;
    e_index[2][1] = a;
  }

  Face& operator=(const Face& that)
  {
    for (int k = 0; k < 3; ++k)
    {
      v_index[k] = that.v_index[k];
      n_index[k] = that.n_index[k];
      vt_index[k] = that.vt_index[k];
      e_index[k][0] = that.e_index[k][0];
      e_index[k][1] = that.e_index[k][1];
      e_index[k][2] = that.e_index[k][2];
      e_index[k][3] = that.e_index[k][3];
    }
    geometric_normal = that.geometric_normal;
    assigned_normal = that.assigned_normal;
    invalidated = that.invalidated;
    return *this;
  }
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
 *  [2] = Partition ID of adjecent cell.
 */
struct PolyFace
{
  std::vector<int> v_indices;
  std::vector<std::vector<int>> edges;
  int face_indices[3];

  Vector3 geometric_normal;
  Vector3 face_centroid;

  bool invalidated;

  PolyFace() { invalidated = false; }
};

} // namespace opensn
