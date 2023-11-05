#include "SurfaceMeshLogicalVolume.h"

#include "mesh/chi_mesh.h"
#include "mesh/SurfaceMesh/chi_surfacemesh.h"
#include "mesh/Raytrace/raytracer.h"

#include "ChiObjectFactory.h"

#include <utility>

namespace chi_mesh
{

RegisterChiObject(chi_mesh, SurfaceMeshLogicalVolume);

chi::InputParameters
SurfaceMeshLogicalVolume::GetInputParameters()
{
  chi::InputParameters params = LogicalVolume::GetInputParameters();

  params.SetDocGroup("LuaLogicVolumes");

  params.AddRequiredParameter<size_t>("surface_mesh_handle",
                                      "Handle to a surface mesh that will represent this object");

  return params;
}

SurfaceMeshLogicalVolume::SurfaceMeshLogicalVolume(const chi::InputParameters& params)
  : LogicalVolume(params),
    surf_mesh(Chi::GetStackItemPtrAsType<chi_mesh::SurfaceMesh>(
      Chi::object_stack, params.GetParamValue<size_t>("surface_mesh_handle"), __FUNCTION__)),
    xbounds_({1.0e6, -1.0e6}),
    ybounds_({1.0e6, -1.0e6}),
    zbounds_({1.0e6, -1.0e6})
{
  const auto& vertices = surf_mesh->GetVertices();
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

// ###################################################################
/**Logical operation for surface mesh.*/
bool
SurfaceMeshLogicalVolume::Inside(const chi_mesh::Vector3& point) const
{
  double tolerance = 1.0e-5;

  //============================================= Boundbox check
  double x = point.x;
  double y = point.y;
  double z = point.z;

  if (not((x >= xbounds_[0]) and (x <= xbounds_[1]))) return false;
  if (not((y >= ybounds_[0]) and (y <= ybounds_[1]))) return false;
  if (not((z >= zbounds_[0]) and (z <= zbounds_[1]))) return false;

  //============================================= Cheapshot pass
  // This pass purely checks if the point have a
  // negative sense with all the faces of the surface.
  // If it does then .. bonus .. we don't need to do
  // anything more because the surface is probably convex.
  bool cheap_pass = true; // now try to disprove
  for (auto& face : surf_mesh->GetTriangles())
  {
    chi_mesh::Vector3 fc = face.face_centroid;
    chi_mesh::Vector3 p_to_fc = fc - point;

    p_to_fc = p_to_fc / p_to_fc.Norm();

    double sense = p_to_fc.Dot(face.geometric_normal);

    if (sense < (0.0 - tolerance))
    {
      cheap_pass = false;
      break;
    }
  } // for f

  // if (!cheap_pass) return false;
  if (cheap_pass) return true;

  //============================================= Expensive pass
  // Getting to here means the cheap pass produced
  // a negative and now we need to do more work.
  for (size_t f = 0; f < surf_mesh->GetTriangles().size(); f++)
  {
    chi_mesh::Vector3 fc = surf_mesh->GetTriangles()[f].face_centroid;
    chi_mesh::Vector3 p_to_fc = fc - point;
    double distance_to_face = p_to_fc.Norm();
    double closest_distance = 1.0e16;
    bool closest_sense_pos = false;

    p_to_fc = p_to_fc / p_to_fc.Norm();

    double sense = p_to_fc.Dot(surf_mesh->GetTriangles()[f].geometric_normal);

    bool good_to_go = true;
    if (sense < (0.0 - tolerance))
    {
      good_to_go = false;
      for (size_t fi = 0; fi < surf_mesh->GetTriangles().size(); fi++)
      {
        if (fi == f) continue; // Skip same face

        // Get all the vertices
        int v0_i = surf_mesh->GetTriangles()[fi].v_index[0];
        int v1_i = surf_mesh->GetTriangles()[fi].v_index[1];
        int v2_i = surf_mesh->GetTriangles()[fi].v_index[2];
        chi_mesh::Vertex v0 = surf_mesh->GetVertices()[v0_i];
        chi_mesh::Vertex v1 = surf_mesh->GetVertices()[v1_i];
        chi_mesh::Vertex v2 = surf_mesh->GetVertices()[v2_i];

        //=========================== Check if the line intersects plane
        chi_mesh::Vertex intp; // Intersection point
        std::pair<double, double> weights;
        bool intersects_plane = chi_mesh::CheckPlaneLineIntersect(
          surf_mesh->GetTriangles()[fi].geometric_normal, v0, point, fc, intp, &weights);
        if (!intersects_plane) continue;

        //=========================== Check if the line intersects the triangle
        bool intersects_triangle = true;

        // Compute the legs
        chi_mesh::Vector3 v01 = v1 - v0;
        chi_mesh::Vector3 v12 = v2 - v1;
        chi_mesh::Vector3 v20 = v0 - v2;

        // Compute the vertices to the point
        chi_mesh::Vector3 v0p = intp - v0;
        chi_mesh::Vector3 v1p = intp - v1;
        chi_mesh::Vector3 v2p = intp - v2;

        // Compute the cross products
        chi_mesh::Vector3 x0p = v01.Cross(v0p);
        chi_mesh::Vector3 x1p = v12.Cross(v1p);
        chi_mesh::Vector3 x2p = v20.Cross(v2p);

        // Normalize them
        x0p = x0p / x0p.Norm();
        x1p = x1p / x1p.Norm();
        x2p = x2p / x2p.Norm();

        chi_mesh::Vector3 face_norm = surf_mesh->GetTriangles()[fi].geometric_normal /
                                      surf_mesh->GetTriangles()[fi].geometric_normal.Norm();

        if (x0p.Dot(face_norm) < 0.0) intersects_triangle = false;
        if (x1p.Dot(face_norm) < 0.0) intersects_triangle = false;
        if (x2p.Dot(face_norm) < 0.0) intersects_triangle = false;

        if (!intersects_triangle) continue;

        //============================ Determine the sense with the triangle
        double sense_with_this_tri = p_to_fc.Dot(surf_mesh->GetTriangles()[fi].geometric_normal);
        double distance_to_triangle = weights.second * distance_to_face;

        if (distance_to_triangle < closest_distance)
        {
          closest_distance = distance_to_triangle;

          if (sense_with_this_tri > 0.0) closest_sense_pos = true;
          else
            closest_sense_pos = false;
        } //

      } // for inner iter face
    }   // if sense negative

    if ((closest_distance < distance_to_face) && closest_sense_pos) good_to_go = true;

    if (!good_to_go) return false;
  } // for f

  return true;
}

} // namespace chi_mesh
