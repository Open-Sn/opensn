// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/logical_volume/surface_mesh_logical_volume.h"
#include "framework/mesh/mesh.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/mesh/raytrace/raytracer.h"
#include "framework/object_factory.h"
#include <utility>

namespace opensn
{

OpenSnRegisterObjectInNamespace(logvol, SurfaceMeshLogicalVolume);

InputParameters
SurfaceMeshLogicalVolume::GetInputParameters()
{
  InputParameters params = LogicalVolume::GetInputParameters();

  params.SetDocGroup("LuaLogicVolumes");

  params.AddRequiredParameter<size_t>("surface_mesh_handle",
                                      "Handle to a surface mesh that will represent this object");

  return params;
}

SurfaceMeshLogicalVolume::SurfaceMeshLogicalVolume(const InputParameters& params)
  : LogicalVolume(params),
    surf_mesh_(GetStackItemPtrAsType<SurfaceMesh>(
      object_stack, params.ParamValue<size_t>("surface_mesh_handle"), __FUNCTION__)),
    xbounds_({1.0e6, -1.0e6}),
    ybounds_({1.0e6, -1.0e6}),
    zbounds_({1.0e6, -1.0e6})
{
  const auto& vertices = surf_mesh_->Vertices();
  for (auto& vertex : vertices)
  {
    const double x = vertex.x;
    const double y = vertex.y;
    const double z = vertex.z;
    if (std::addressof(vertex) == std::addressof(vertices.front()))
    {
      xbounds_[0] = x;
      xbounds_[1] = x;
      ybounds_[0] = y;
      ybounds_[1] = y;
      zbounds_[0] = z;
      zbounds_[1] = z;
    }
    else
    {
      xbounds_[0] = std::min(xbounds_[0], x);
      xbounds_[1] = std::max(xbounds_[1], x);
      ybounds_[0] = std::min(ybounds_[0], y);
      ybounds_[1] = std::max(ybounds_[1], y);
      zbounds_[0] = std::min(zbounds_[0], z);
      zbounds_[1] = std::max(zbounds_[1], z);
    }
  }
}

bool
SurfaceMeshLogicalVolume::Inside(const Vector3& point) const
{
  double tolerance = 1.0e-5;

  // Boundbox check
  double x = point.x;
  double y = point.y;
  double z = point.z;

  if (not((x >= xbounds_[0]) and (x <= xbounds_[1])))
    return false;
  if (not((y >= ybounds_[0]) and (y <= ybounds_[1])))
    return false;
  if (not((z >= zbounds_[0]) and (z <= zbounds_[1])))
    return false;

  // Cheapshot pass
  // This pass purely checks if the point have a
  // negative sense with all the faces of the surface.
  // If it does then .. bonus .. we don't need to do
  // anything more because the surface is probably convex.
  bool cheap_pass = true; // now try to disprove
  for (auto& face : surf_mesh_->Triangles())
  {
    Vector3 fc = face.face_centroid;
    Vector3 p_to_fc = fc - point;

    p_to_fc = p_to_fc / p_to_fc.Norm();

    double sense = p_to_fc.Dot(face.geometric_normal);

    if (sense < (0.0 - tolerance))
    {
      cheap_pass = false;
      break;
    }
  } // for f

  // if (!cheap_pass) return false;
  if (cheap_pass)
    return true;

  // Expensive pass
  // Getting to here means the cheap pass produced
  // a negative and now we need to do more work.
  for (size_t f = 0; f < surf_mesh_->Triangles().size(); ++f)
  {
    Vector3 fc = surf_mesh_->Triangles()[f].face_centroid;
    Vector3 p_to_fc = fc - point;
    double distance_to_face = p_to_fc.Norm();
    double closest_distance = 1.0e16;
    bool closest_sense_pos = false;

    p_to_fc = p_to_fc / p_to_fc.Norm();

    double sense = p_to_fc.Dot(surf_mesh_->Triangles()[f].geometric_normal);

    bool good_to_go = true;
    if (sense < (0.0 - tolerance))
    {
      good_to_go = false;
      for (size_t fi = 0; fi < surf_mesh_->Triangles().size(); ++fi)
      {
        if (fi == f)
          continue; // Skip same face

        // Get all the vertices
        int v0_i = surf_mesh_->Triangles()[fi].v_index[0];
        int v1_i = surf_mesh_->Triangles()[fi].v_index[1];
        int v2_i = surf_mesh_->Triangles()[fi].v_index[2];
        Vector3 v0 = surf_mesh_->Vertices()[v0_i];
        Vector3 v1 = surf_mesh_->Vertices()[v1_i];
        Vector3 v2 = surf_mesh_->Vertices()[v2_i];

        // Check if the line intersects plane
        Vector3 intp; // Intersection point
        std::pair<double, double> weights;
        bool intersects_plane = CheckPlaneLineIntersect(
          surf_mesh_->Triangles()[fi].geometric_normal, v0, point, fc, intp, &weights);
        if (not intersects_plane)
          continue;

        // Check if the line intersects the triangle
        bool intersects_triangle = true;

        // Compute the legs
        Vector3 v01 = v1 - v0;
        Vector3 v12 = v2 - v1;
        Vector3 v20 = v0 - v2;

        // Compute the vertices to the point
        Vector3 v0p = intp - v0;
        Vector3 v1p = intp - v1;
        Vector3 v2p = intp - v2;

        // Compute the cross products
        Vector3 x0p = v01.Cross(v0p);
        Vector3 x1p = v12.Cross(v1p);
        Vector3 x2p = v20.Cross(v2p);

        // Normalize them
        x0p = x0p / x0p.Norm();
        x1p = x1p / x1p.Norm();
        x2p = x2p / x2p.Norm();

        Vector3 face_norm = surf_mesh_->Triangles()[fi].geometric_normal /
                            surf_mesh_->Triangles()[fi].geometric_normal.Norm();

        if (x0p.Dot(face_norm) < 0.0)
          intersects_triangle = false;
        if (x1p.Dot(face_norm) < 0.0)
          intersects_triangle = false;
        if (x2p.Dot(face_norm) < 0.0)
          intersects_triangle = false;

        if (not intersects_triangle)
          continue;

        // Determine the sense with the triangle
        double sense_with_this_tri = p_to_fc.Dot(surf_mesh_->Triangles()[fi].geometric_normal);
        double distance_to_triangle = weights.second * distance_to_face;

        if (distance_to_triangle < closest_distance)
        {
          closest_distance = distance_to_triangle;

          if (sense_with_this_tri > 0.0)
            closest_sense_pos = true;
          else
            closest_sense_pos = false;
        } //

      } // for inner iter face
    }   // if sense negative

    if ((closest_distance < distance_to_face) and closest_sense_pos)
      good_to_go = true;

    if (not good_to_go)
      return false;
  } // for f

  return true;
}

} // namespace opensn
