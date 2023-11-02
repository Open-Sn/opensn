#include "chi_ffinter_slice.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_mpi.h"
#include "chi_log.h"

#include <fstream>

/***/
void chi_mesh::FieldFunctionInterpolationSlice::ExportPython(std::string base_name)
{
  std::ofstream ofile;

  std::string fileName = base_name;
  fileName = fileName + std::to_string(Chi::mpi.location_id);
  fileName = fileName + std::string(".py");
  ofile.open(fileName);

  ofile
  << "import numpy as np\n"
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

  ofile
  << "class CellData:\n"
     "  def __init__(self):\n"
     "    self.data_points=[]\n"
     "    self.xy = []\n"
     "    self.c = []\n\n";

  std::string offset;
  if (Chi::mpi.location_id == 0)
  {
    std::string submod_name = base_name;
    submod_name = submod_name + std::to_string(Chi::mpi.location_id+1);

    if (Chi::mpi.process_count>1)
    {
      ofile << "import " << submod_name << "\n\n";
    }

    ofile
    << "class BaseDataClass:\n"
    << "  def __init__(self):\n"
    << "    data_object = []\n"
    << "    self.data_object = data_object\n";

    offset = std::string("    ");
  }
  else if (Chi::mpi.process_count>1)
  {

    if (Chi::mpi.location_id != (Chi::mpi.process_count-1))
    {
      std::string submod_name = base_name;
      submod_name = submod_name + std::to_string(Chi::mpi.location_id+1);

      ofile << "import " << submod_name << "\n\n";
    }

    ofile
      << "def AddData(data_object):\n";

    offset = std::string("  ");
  }

  size_t num_cells = cell_intersections_.size();
  for (int c=0; c<num_cells; c++)
  {
    double x = 0.0;
    double y = 0.0;
    double v = 0.0;

    ofile
    << offset << "new_cell_data = CellData()\n";

    size_t num_points = cell_intersections_[c].intersections.size();
    ofile
    << offset << "new_cell_data.xy = np.zeros(("
              << std::to_string(num_points) << ",2))\n"
    << offset << "new_cell_data.c = np.zeros("
              << std::to_string(num_points) << ")\n";
    for (int p=0; p<num_points; p++)
    {
      x = cell_intersections_[c].intersections[p].point2d.x;
      y = cell_intersections_[c].intersections[p].point2d.y;
      v = cell_intersections_[c].intersections[p].point_value;

      ofile
        << offset
        << "new_cell_data.xy[" << std::to_string(p)
        << ",0] = " << std::to_string(x) << "\n"
        << offset
        << "new_cell_data.xy[" << std::to_string(p)
        << ",1] = " << std::to_string(y) << "\n"
        << offset
        << "new_cell_data.c[" << std::to_string(p)
        << "] = " << std::to_string(v) << "\n";
    }
    v = cell_intersections_[c].cell_avg_value;
    ofile
      << offset
      << "new_cell_data.avg = " << std::to_string(v) << "\n";


    ofile
    << offset << "data_object.append(new_cell_data)\n";
  }

  if (Chi::mpi.location_id != (Chi::mpi.process_count-1))
  {
    std::string submod_name = base_name;
    submod_name = submod_name + std::to_string(Chi::mpi.location_id+1);

    ofile
    << offset << "data_object = "
    << submod_name << ".AddData(data_object)\n\n";
  }
  if (Chi::mpi.location_id>0)
  {
    ofile
    << offset << "return data_object\n";
  }

  if (Chi::mpi.location_id==0)
  {
    ofile
    << "data = BaseDataClass()\n"
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

  Chi::log.Log()
    << "Exported Python files for field func \""
    << field_functions_[0]->TextName()
    << "\" to base name \""
    << base_name << "\" Successfully";


}