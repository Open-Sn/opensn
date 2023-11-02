#include "chi_ffinter_slice.h"
#include "mesh/Cell/cell.h"

#include "mesh/Raytrace/raytracing.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

/**Initializes the data structures necessary for interpolation. This is
 * independent of the physics and hence is a routine on its own.
 *
 * The first step of this initialization is to determine which cells
 * are intersected by this plane. For polyhedrons this is evaluated
 * tet-by-tet.
 *
 * The second step is find where face-edges are intersected. This will
 * effectively create intersection polygons.*/
void
chi_mesh::FieldFunctionInterpolationSlice::Initialize()
{
  Chi::log.Log0Verbose1() << "Initializing slice interpolator.";
  //================================================== Check grid available
  if (field_functions_.empty())
    throw std::logic_error("Unassigned field function in slice "
                           "field function interpolator.");

  const auto& grid = field_functions_.front()->GetSpatialDiscretization().Grid();

  //================================================== Find cells intersecting
  //                                                   plane
  std::vector<uint64_t> intersecting_cell_indices;

  for (const auto& cell : grid.local_cells)
  {
    auto cell_local_index = cell.local_id_;

    if (cell.Type() == chi_mesh::CellType::SLAB)
      throw std::logic_error("FieldFunctionInterpolationSlice "
                             "does not support 1D cells.");
    if (cell.Type() == chi_mesh::CellType::POLYGON)
      intersecting_cell_indices.push_back(cell_local_index);
    else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      bool intersects = false;

      size_t num_faces = cell.faces_.size();
      for (size_t f = 0; f < num_faces; f++)
      {
        const auto& face = cell.faces_[f];

        size_t num_edges = face.vertex_ids_.size();
        for (size_t e = 0; e < num_edges; e++)
        {
          size_t ep1 = (e < (num_edges - 1)) ? e + 1 : 0;
          uint64_t v0_i = face.vertex_ids_[e];
          uint64_t v1_i = face.vertex_ids_[ep1];

          std::vector<chi_mesh::Vector3> tet_points;

          tet_points.push_back(grid.vertices[v0_i]);
          tet_points.push_back(grid.vertices[v1_i]);
          tet_points.push_back(cell.faces_[f].centroid_);
          tet_points.push_back(cell.centroid_);

          if (CheckPlaneTetIntersect(this->normal_, this->plane_point_, tet_points))
          {
            intersecting_cell_indices.push_back(cell_local_index);
            intersects = true;
            break; // from for e
          }
        }                      // for e
        if (intersects) break; // from for f
      }                        // for f
    }
    else
      throw std::logic_error("Unsupported cell type in call "
                             "to Slice Initialize.");
  } // for local cell

  //================================================== Computing cell
  //                                                   intersections
  for (const uint64_t cell_local_index : intersecting_cell_indices)
  {
    const auto& cell = grid.local_cells[cell_local_index];

    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      FFICellIntersection cell_isds;
      cell_isds.ref_cell_local_id = cell_local_index;

      //========================================= Loop over vertices
      for (uint64_t v0gi : cell.vertex_ids_)
      {
        FFIFaceEdgeIntersection face_isds;

        const auto nudge = 1.0e-4 * (grid.vertices[v0gi] - cell.centroid_);

        face_isds.point = grid.vertices[v0gi] - nudge;
        face_isds.point2d = grid.vertices[v0gi];
        cell_isds.intersections.push_back(face_isds);
      }

      //========================================= Set intersection center
      cell_isds.intersection_centre = cell.centroid_;

      //========================================= Set straight 2D center
      // This is normally transformed for the 3D case
      cell_isds.intersection_2d_centre = cell.centroid_;

      //========================================= Same for 2D points
      size_t num_points = cell_isds.intersections.size();
      for (size_t p = 0; p < num_points; p++)
      {
        chi_mesh::Vector3 vref = cell_isds.intersections[p].point - plane_point_;

        cell_isds.intersections[p].point2d = vref;
      }

      cell_intersections_.push_back(cell_isds);
    } // polygon
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      //========================================= Initialize cell intersection
      //                                          data structure
      FFICellIntersection cell_isds;
      cell_isds.ref_cell_local_id = cell_local_index;

      //========================================= Loop over faces
      size_t num_faces = cell.faces_.size();
      for (size_t f = 0; f < num_faces; f++)
      {
        const auto& face = cell.faces_[f];
        //================================== Loop over edges
        size_t num_edges = face.vertex_ids_.size();
        for (size_t e = 0; e < num_edges; e++)
        {
          size_t ep1 = (e < (num_edges - 1)) ? e + 1 : 0;
          uint64_t v0gi = face.vertex_ids_[e];   // global index v0
          uint64_t v1gi = face.vertex_ids_[ep1]; // global index v1

          const auto& v0 = grid.vertices[v0gi];
          const auto& v1 = grid.vertices[v1gi];

          chi_mesh::Vertex interstion_point; // Placeholder
          std::pair<double, double> weights;

          //=========================== Check if intersects plane
          if (CheckPlaneLineIntersect(
                this->normal_, this->plane_point_, v0, v1, interstion_point, &weights))
          {
            //==================== Check for duplicate
            bool duplicate_found = false;
            for (auto& existing_face_isds : cell_isds.intersections)
            {
              double dif = (existing_face_isds.point - interstion_point).NormSquare();
              if (dif < 1.0e-6)
              {
                duplicate_found = true;
                break;
              }
            }

            //==================== No duplicate
            if (!duplicate_found)
            {
              FFIFaceEdgeIntersection face_isds;

              face_isds.point = interstion_point;
              cell_isds.intersections.push_back(face_isds);
            }
          } // if intersecting
        }   // for edge
      }     // for face

      //====================================

      //==================================== Computing intersection centre
      size_t num_points = cell_isds.intersections.size();
      if (num_points > 0)
      {
        for (int p = 0; p < num_points; p++)
        {
          cell_isds.intersection_centre =
            cell_isds.intersection_centre + cell_isds.intersections[p].point;
        }
        cell_isds.intersection_centre =
          cell_isds.intersection_centre / static_cast<double>(num_points);
      }
      else
      {
        Chi::log.LogAllWarning() << "No face intersections encountered "
                                    "for a cell that is indicated as being "
                                    "intersected. Slice FF interp.";
      }

      //==================================== Computing 2D transforms
      chi_mesh::Vector3 vref = cell_isds.intersection_centre - plane_point_;

      cell_isds.intersection_2d_centre.x = vref.Dot(tangent_);
      cell_isds.intersection_2d_centre.y = vref.Dot(binorm_);
      cell_isds.intersection_2d_centre.z = vref.Dot(normal_);

      //==================================== Points
      std::vector<FFIFaceEdgeIntersection> unsorted_points;
      for (int p = 0; p < num_points; p++)
      {
        vref = cell_isds.intersections[p].point - plane_point_;

        cell_isds.intersections[p].point2d.x = vref.Dot(tangent_);
        cell_isds.intersections[p].point2d.y = vref.Dot(binorm_);
        cell_isds.intersections[p].point2d.z = vref.Dot(normal_);

        unsorted_points.push_back(cell_isds.intersections[p]);
      }
      cell_isds.intersections.clear();

      //==================================== Sort points clockwise
      // The first point is retrieved from the unused stack.
      // Subsequent points are only added if they form a
      // convex line wrt the right hand rule.
      cell_isds.intersections.push_back(unsorted_points[0]);
      unsorted_points.erase(unsorted_points.begin());

      while (!unsorted_points.empty())
      {
        for (int p = 0; p < unsorted_points.size(); p++)
        {
          chi_mesh::Vector3 v1 =
            unsorted_points[p].point2d - cell_isds.intersections.back().point2d;

          bool illegal_value = false;
          for (int pr = 0; pr < unsorted_points.size(); pr++)
          {
            if (pr != p)
            {
              chi_mesh::Vector3 vr = unsorted_points[pr].point2d - unsorted_points[p].point2d;

              if (vr.Cross(v1).z < 0.0)
              {
                illegal_value = true;
                break;
              }
            } // if not p
          }   // for pr

          if (!illegal_value)
          {
            cell_isds.intersections.push_back(unsorted_points[p]);
            unsorted_points.erase(unsorted_points.begin() + p);
            break;
          }

        } // for p
      }
      cell_intersections_.push_back(cell_isds);
    } // polyhedron
  }   // for intersected cell

  // chi::log.Log() << "Finished initializing interpolator.";
}