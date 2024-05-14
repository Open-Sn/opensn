// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/interpolation/ffinter_slice.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/raytrace/raytracer.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/vector_ghost_communicator/vector_ghost_communicator.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <fstream>

namespace opensn
{

void
FieldFunctionInterpolationSlice::Initialize()
{
  log.Log0Verbose1() << "Initializing slice interpolator.";
  // Check grid available
  if (field_functions_.empty())
    throw std::logic_error("Unassigned field function in slice "
                           "field function interpolator.");

  const auto& grid = field_functions_.front()->GetSpatialDiscretization().Grid();

  // Find cells intersecting plane
  std::vector<uint64_t> intersecting_cell_indices;

  for (const auto& cell : grid.local_cells)
  {
    auto cell_local_index = cell.local_id_;

    if (cell.Type() == CellType::SLAB)
      throw std::logic_error("FieldFunctionInterpolationSlice "
                             "does not support 1D cells.");
    if (cell.Type() == CellType::POLYGON)
      intersecting_cell_indices.push_back(cell_local_index);
    else if (cell.Type() == CellType::POLYHEDRON)
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

          std::vector<Vector3> tet_points;

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
        } // for e
        if (intersects)
          break; // from for f
      }          // for f
    }
    else
      throw std::logic_error("Unsupported cell type in call "
                             "to Slice Initialize.");
  } // for local cell

  // Computing cell intersections
  for (const uint64_t cell_local_index : intersecting_cell_indices)
  {
    const auto& cell = grid.local_cells[cell_local_index];

    if (cell.Type() == CellType::POLYGON)
    {
      FFICellIntersection cell_isds;
      cell_isds.ref_cell_local_id = cell_local_index;

      // Loop over vertices
      for (uint64_t v0gi : cell.vertex_ids_)
      {
        FFIFaceEdgeIntersection face_isds;

        const auto nudge = 1.0e-4 * (grid.vertices[v0gi] - cell.centroid_);

        face_isds.point = grid.vertices[v0gi] - nudge;
        face_isds.point2d = grid.vertices[v0gi];
        cell_isds.intersections.push_back(face_isds);
      }

      // Set intersection center
      cell_isds.intersection_centre = cell.centroid_;

      // Set straight 2D center
      // This is normally transformed for the 3D case
      cell_isds.intersection_2d_centre = cell.centroid_;

      // Same for 2D points
      size_t num_points = cell_isds.intersections.size();
      for (size_t p = 0; p < num_points; p++)
      {
        Vector3 vref = cell_isds.intersections[p].point - plane_point_;

        cell_isds.intersections[p].point2d = vref;
      }

      cell_intersections_.push_back(cell_isds);
    } // polygon
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    else if (cell.Type() == CellType::POLYHEDRON)
    {
      // Initialize cell intersection data structure
      FFICellIntersection cell_isds;
      cell_isds.ref_cell_local_id = cell_local_index;

      // Loop over faces
      size_t num_faces = cell.faces_.size();
      for (size_t f = 0; f < num_faces; f++)
      {
        const auto& face = cell.faces_[f];
        // Loop over edges
        size_t num_edges = face.vertex_ids_.size();
        for (size_t e = 0; e < num_edges; e++)
        {
          size_t ep1 = (e < (num_edges - 1)) ? e + 1 : 0;
          uint64_t v0gi = face.vertex_ids_[e];   // global index v0
          uint64_t v1gi = face.vertex_ids_[ep1]; // global index v1

          const auto& v0 = grid.vertices[v0gi];
          const auto& v1 = grid.vertices[v1gi];

          Vertex interstion_point; // Placeholder
          std::pair<double, double> weights;

          // Check if intersects plane
          if (CheckPlaneLineIntersect(
                this->normal_, this->plane_point_, v0, v1, interstion_point, &weights))
          {
            // Check for duplicate
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

            // No duplicate
            if (not duplicate_found)
            {
              FFIFaceEdgeIntersection face_isds;

              face_isds.point = interstion_point;
              cell_isds.intersections.push_back(face_isds);
            }
          } // if intersecting
        }   // for edge
      }     // for face

      // Computing intersection centre
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
        log.LogAllWarning() << "No face intersections encountered "
                               "for a cell that is indicated as being "
                               "intersected. Slice FF interp.";
      }

      // Computing 2D transforms
      Vector3 vref = cell_isds.intersection_centre - plane_point_;

      cell_isds.intersection_2d_centre.x = vref.Dot(tangent_);
      cell_isds.intersection_2d_centre.y = vref.Dot(binorm_);
      cell_isds.intersection_2d_centre.z = vref.Dot(normal_);

      // Points
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

      // Sort points clockwise
      // The first point is retrieved from the unused stack.
      // Subsequent points are only added if they form a
      // convex line wrt the right hand rule.
      cell_isds.intersections.push_back(unsorted_points[0]);
      unsorted_points.erase(unsorted_points.begin());

      while (not unsorted_points.empty())
      {
        for (int p = 0; p < unsorted_points.size(); p++)
        {
          Vector3 v1 = unsorted_points[p].point2d - cell_isds.intersections.back().point2d;

          bool illegal_value = false;
          for (int pr = 0; pr < unsorted_points.size(); pr++)
          {
            if (pr != p)
            {
              Vector3 vr = unsorted_points[pr].point2d - unsorted_points[p].point2d;

              if (vr.Cross(v1).z < 0.0)
              {
                illegal_value = true;
                break;
              }
            } // if not p
          }   // for pr

          if (not illegal_value)
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

void
FieldFunctionInterpolationSlice::Execute()
{
  const auto& ref_ff = *field_functions_.front();
  const auto& sdm = ref_ff.GetSpatialDiscretization();
  const auto& grid = sdm.Grid();

  const auto& uk_man = ref_ff.GetUnknownManager();
  const auto uid = 0;
  const auto cid = ref_component_;

  const auto field_data = ref_ff.GetGhostedFieldVector();

  for (auto& cell_intersection : cell_intersections_)
  {
    const auto& cell = grid.local_cells[cell_intersection.ref_cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    std::vector<double> dof_values(num_nodes, 0.0);
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);
      dof_values[i] = field_data[imap];
    }

    std::vector<double> shape_values(num_nodes, 0.0);
    for (auto& edge_intersection : cell_intersection.intersections)
    {
      cell_mapping.ShapeValues(edge_intersection.point, shape_values);
      double point_value = 0.0;
      for (size_t i = 0; i < num_nodes; ++i)
        point_value += dof_values[i] * shape_values[i];

      edge_intersection.point_value = point_value;
    } // for edge intersection
  }   // for cell intersection
}

void
FieldFunctionInterpolationSlice::ExportToPython(std::string base_name)
{
  std::ofstream ofile;

  std::string fileName = base_name;
  fileName = fileName + std::to_string(opensn::mpi_comm.rank());
  fileName = fileName + std::string(".py");
  ofile.open(fileName);

  ofile << "import numpy as np\n"
           "import matplotlib.pyplot as plt\n"
           "import matplotlib.cm as cm\n"
           "from matplotlib.collections import PatchCollection\n"
           "from matplotlib.patches import Polygon\n"
           "import matplotlib.colors as colors\n"
        << "\n"
        << "class Datapoint:\n"
           "  def __init__(self,x,y,value):\n"
           "    self.x = x\n"
           "    self.y = y\n"
           "    self.value = value\n\n";

  ofile << "class CellData:\n"
           "  def __init__(self):\n"
           "    self.data_points=[]\n"
           "    self.xy = []\n"
           "    self.c = []\n\n";

  std::string offset;
  if (opensn::mpi_comm.rank() == 0)
  {
    std::string submod_name = base_name;
    submod_name = submod_name + std::to_string(opensn::mpi_comm.rank() + 1);

    if (opensn::mpi_comm.size() > 1)
    {
      ofile << "import " << submod_name << "\n\n";
    }

    ofile << "class BaseDataClass:\n"
          << "  def __init__(self):\n"
          << "    data_object = []\n"
          << "    self.data_object = data_object\n";

    offset = std::string("    ");
  }
  else if (opensn::mpi_comm.size() > 1)
  {

    if (opensn::mpi_comm.rank() != (opensn::mpi_comm.size() - 1))
    {
      std::string submod_name = base_name;
      submod_name = submod_name + std::to_string(opensn::mpi_comm.rank() + 1);

      ofile << "import " << submod_name << "\n\n";
    }

    ofile << "def AddData(data_object):\n";

    offset = std::string("  ");
  }

  size_t num_cells = cell_intersections_.size();
  for (int c = 0; c < num_cells; c++)
  {
    double x = 0.0;
    double y = 0.0;
    double v = 0.0;

    ofile << offset << "new_cell_data = CellData()\n";

    size_t num_points = cell_intersections_[c].intersections.size();
    ofile << offset << "new_cell_data.xy = np.zeros((" << std::to_string(num_points) << ",2))\n"
          << offset << "new_cell_data.c = np.zeros(" << std::to_string(num_points) << ")\n";
    for (int p = 0; p < num_points; p++)
    {
      x = cell_intersections_[c].intersections[p].point2d.x;
      y = cell_intersections_[c].intersections[p].point2d.y;
      v = cell_intersections_[c].intersections[p].point_value;

      ofile << offset << "new_cell_data.xy[" << std::to_string(p) << ",0] = " << std::to_string(x)
            << "\n"
            << offset << "new_cell_data.xy[" << std::to_string(p) << ",1] = " << std::to_string(y)
            << "\n"
            << offset << "new_cell_data.c[" << std::to_string(p) << "] = " << std::to_string(v)
            << "\n";
    }
    v = cell_intersections_[c].cell_avg_value;
    ofile << offset << "new_cell_data.avg = " << std::to_string(v) << "\n";

    ofile << offset << "data_object.append(new_cell_data)\n";
  }

  if (opensn::mpi_comm.rank() != (opensn::mpi_comm.size() - 1))
  {
    std::string submod_name = base_name;
    submod_name = submod_name + std::to_string(opensn::mpi_comm.rank() + 1);

    ofile << offset << "data_object = " << submod_name << ".AddData(data_object)\n\n";
  }
  if (opensn::mpi_comm.rank() > 0)
  {
    ofile << offset << "return data_object\n";
  }

  if (opensn::mpi_comm.rank() == 0)
  {
    ofile << "data = BaseDataClass()\n"
          << "print(len(data.data_object))\n\n"
          << "N = len(data.data_object)\n"
             "\n"
             "maxavg = 0.0\n"
             "maxc = 0.0\n"
             "xycount = 0\n"
             "for c in range(0,N):\n"
             "    for i in range(0,np.size(data.data_object[c].c)):\n"
             "        xycount += 1\n"
             "        if data.data_object[c].avg>maxavg:\n"
             "            maxavg = data.data_object[c].avg\n"
             "        if data.data_object[c].c[i]>maxc:\n"
             "            maxc = data.data_object[c].c[i]\n"
             "            \n"
             "xmax = -9999999.0\n"
             "xmin =  9999999.0\n"
             "ymax = -9999999.0\n"
             "ymin =  9999999.0\n"
             "zmax = -9999999.0\n"
             "zmin =  9999999.0\n"
             "x = np.zeros(xycount)\n"
             "y = np.zeros(xycount)\n"
             "z = np.zeros(xycount)\n"
             "xycount = -1\n"
             "for c in range(0,N):\n"
             "    for i in range(0,np.size(data.data_object[c].c)):\n"
             "        xycount += 1\n"
             "        x[xycount] = data.data_object[c].xy[i,0]\n"
             "        y[xycount] = data.data_object[c].xy[i,1]\n"
             "        z[xycount] = data.data_object[c].c[i]\n"
             "\n"
             "        if x[xycount]>xmax: xmax = x[xycount]\n"
             "        if x[xycount]<xmin: xmin = x[xycount]\n"
             "        if y[xycount]>ymax: ymax = y[xycount]\n"
             "        if y[xycount]<ymin: ymin = y[xycount]\n"
             "        if z[xycount]>zmax: zmax = z[xycount]\n"
             "        if z[xycount]<zmin: zmin = z[xycount]\n"
             "\n"
             "print(\"phi_max=%g phi_min=%g\" %(zmax,zmin))\n"
             "\n"
             "fig,ax = plt.subplots(1)\n"
             "cmapjet = plt.get_cmap('jet')\n"
             "cNorm = colors.Normalize(vmin=0,vmax=maxavg)\n"
             "scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmapjet)\n"
             "\n"
             "cb_scale = np.linspace(zmin,zmax*1.00001, 124, endpoint=True)\n"
             "\n"
             "cntr1 = plt.tricontourf(x,y,z,cb_scale,cmap=cmapjet)\n"
             "\n"
             "for c in range(0,N):\n"
             "    col = scalarMap.to_rgba(data.data_object[c].avg)\n"
             "    poly = Polygon(data.data_object[c].xy,\n"
             "                   closed=True,linestyle='-',fill=False)\n"
             "    patch = []\n"
             "    patch.append(poly)\n"
             "    coll = PatchCollection(patch)\n"
             "    coll.set_facecolor([0,0,0,0])\n"
             "    coll.set_edgecolor([0,0,0,1])\n"
             "    coll.set_linewidth(0.3)\n"
             "\n"
             "    ax.add_collection(coll)\n"
             "\n"
             "cb = fig.colorbar(cntr1,ax=ax)\n"
             "cb.set_ticks(np.linspace(zmin,zmax, 11, endpoint=True))\n"
             "ax.set_xlim([xmin,xmax])\n"
             "ax.set_ylim([ymin,ymax])\n"
             "plt.show()\n";
  }

  ofile.close();

  log.Log() << "Exported Python files for field func \"" << field_functions_[0]->TextName()
            << "\" to base name \"" << base_name << "\" Successfully";
}

} // namespace opensn
