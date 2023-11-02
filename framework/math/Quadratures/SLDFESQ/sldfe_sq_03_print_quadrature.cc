#include "sldfe_sq.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Prints the quadrature to file.*/
void
chi_math::SimplifiedLDFESQ::Quadrature::PrintQuadratureToFile()
{
  Chi::log.Log() << "Printing SLDFE-Quadrature to file.";

  std::ofstream vert_file, cell_file, points_file, python_file;
  vert_file.open(output_filename_prefix_ + "verts.txt");
  {
    for (const auto& sq : deployed_SQs_)
      for (int v = 0; v < 4; ++v)
      {
        auto& v0 = sq.vertices_xyz_prime[v];
        auto& v1 = sq.vertices_xyz_prime[(v < 3) ? v + 1 : 0];

        for (int d = 0; d <= 10; ++d)
        {
          auto vert = (1.0 - d / 10.0) * v0 + (d / 10.0) * v1;
          vert = vert * sq.octant_modifier;
          vert.Normalize();
          vert_file << vert.x << " " << vert.y << " " << vert.z << "\n";
        }
      }
  }
  vert_file.close();

  cell_file.open(output_filename_prefix_ + "cells.txt");
  {
    int vi = 0;
    for (const auto& sq : deployed_SQs_)
    {
      for (const auto& vert : sq.vertices_xyz)
        for (int d = 0; d <= 10; ++d)
          cell_file << vi++ << " ";
      cell_file << "\n";
    }
  }
  cell_file.close();

  double total_weight = 0.0;
  points_file.open(output_filename_prefix_ + "points.txt");
  {
    for (auto& sq : deployed_SQs_)
    {
      int ss = -1;
      for (const auto& point : sq.sub_sqr_points)
      {
        ++ss;
        for (int i = 0; i < 3; ++i)
          points_file << point[i] << " ";
        points_file << sq.sub_sqr_weights[ss];
        total_weight += sq.sub_sqr_weights[ss];
        points_file << "\n";
      }
    }
  }
  points_file.close();

  python_file.open(output_filename_prefix_ + "python.py");
  python_file << "import matplotlib.pyplot as plt\n"
                 "from mpl_toolkits import mplot3d\n"
                 "import mpl_toolkits.mplot3d.art3d as art3d\n"
                 "import mpl_toolkits.mplot3d as ax3\n"
                 "import matplotlib.transforms as mpltransform\n"
                 "\n"
                 "import numpy as np\n"
                 "import math\n"
                 "\n"
                 "#====================================== Read vertices\n"
                 "verts = []\n"
                 "verts_file = open(\""
              << output_filename_prefix_
              << "verts.txt\")\n"
                 "for line in verts_file:\n"
                 "    words = line.split()\n"
                 "    "
                 "verts.append(np.array([float(words[0]),float(words[1]),float("
                 "words[2])]))\n"
                 "verts_file.close()\n"
                 "\n"
                 "#====================================== Read cells\n"
                 "cells = []\n"
                 "cells_file = open(\""
              << output_filename_prefix_
              << "cells.txt\")\n"
                 "for line in cells_file:\n"
                 "    words = line.split()\n"
                 "    cell = []\n"
                 "    for word in words:\n"
                 "        cell.append(int(word))\n"
                 "    cells.append(cell)\n"
                 "cells_file.close()\n"
                 "\n"
                 "#====================================== Read points\n"
                 "points = []\n"
                 "weightsum=0.0\n"
                 "points_file = open(\""
              << output_filename_prefix_
              << "points.txt\")\n"
                 "for line in points_file:\n"
                 "    words = line.split()\n"
                 "    point = []\n"
                 "    for word in words:\n"
                 "        point.append(float(word))\n"
                 "    points.append(point)\n"
                 "    weightsum += point[3]\n"
                 "points_file.close()\n"
                 "\n"
                 "print(\"Weightsum check: \",weightsum,weightsum/4/math.pi)\n"
                 "\n"
                 "points_array = np.array(points)\n"
                 "\n"
                 "#====================================== Generate polygons\n"
                 "patches = []\n"
                 "for cell in cells:\n"
                 "\n"
                 "    vertex_list = []\n"
                 "    for index in cell:\n"
                 "        vertex_list.append(verts[index])\n"
                 "\n"
                 "    polygon = art3d.Poly3DCollection([vertex_list])\n"
                 "    polygon.set_color([1.0,1.0,1.0,1.0])\n"
                 "    polygon.set_edgecolor([0.0,0.0,0.0,1.0])\n"
                 "    patches.append(polygon)\n"
                 "\n"
                 "#====================================== Plot polygons\n"
                 "fig = plt.figure(figsize=(10,8.5))\n"
                 "ax = ax3.Axes3D(fig, proj_type = 'ortho')\n"
                 "\n"
                 "ax.view_init(20,25)\n"
                 "limit = 1\n"
                 "end = int(len(patches)/limit)\n"
                 "for poly in patches[0:end]:\n"
                 "  ax.add_collection3d(poly)\n"
                 "\n"
                 "avg_weight = 0.5*math.pi/len(points)\n"
                 "\n"
                 "psize=min(160.0,160*(1.0/avg_weight)*(48/len(points)))\n"
                 "# psize=160\n"
                 "# print(len(points))\n"
                 "end = int(len(points_array)/limit)\n"
                 "# ax.scatter3D(points_array[0:end,0],\n"
                 "#              points_array[0:end,1],\n"
                 "#              points_array[0:end,2],depthshade=False,\n"
                 "#              s=psize*points_array[:,3],c=[[0,0,0,1]])\n"
                 "\n"
                 "\n"
                 "if limit==8:\n"
                 "    ax.set_xlim([0.0,1.0])\n"
                 "    ax.set_ylim([0.0,1.0])\n"
                 "    ax.set_zlim([0.0,1.0])\n"
                 "else:\n"
                 "    ax.set_xlim([-1.0,1.0])\n"
                 "    ax.set_ylim([-1.0,1.0])\n"
                 "    ax.set_zlim([-1.0,1.0])\n"
                 "\n"
                 "ax.margins(0.5)\n"
                 "ax.set_xlabel(r\"$\\mu$\")\n"
                 "ax.set_ylabel(r\"$\\eta$\")\n"
                 "ax.set_zlabel(r\"$\\xi$\")\n"
                 "plt.show()\n";
  python_file.close();

  Chi::log.Log() << "Done printing SLDFE-Quadrature to file.";
}