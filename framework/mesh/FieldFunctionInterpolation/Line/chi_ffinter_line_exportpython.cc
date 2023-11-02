#include "chi_ffinter_line.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_mpi.h"
#include "chi_log.h"

#include <fstream>

//###################################################################
/***/
void chi_mesh::FieldFunctionInterpolationLine::
ExportPython(std::string base_name)
{
  std::ofstream ofile;

  std::string fileName = base_name;
  fileName = fileName + std::to_string(Chi::mpi.location_id);
  fileName = fileName + std::string(".py");
  ofile.open(fileName);

  ofile
    << "import numpy as np\n"
       "import matplotlib.pyplot as plt\n"
    << "\n";

  std::string offset;
  std::string submod_name;
  if (Chi::mpi.location_id == 0)
  {
    submod_name = base_name;
    submod_name = submod_name + std::to_string(Chi::mpi.location_id+1);

    if (Chi::mpi.process_count>1)
    {
      ofile << "import " << submod_name << "\n\n";
    }

    for (int ff=0; ff < field_functions_.size(); ff++)
    {
      ofile
        << "data" << ff << "=np.zeros([" << interpolation_points_.size()
        << ",5])\n";
    }
    for (int ca=0; ca < custom_arrays_.size(); ca++)
    {
      int ff = ca + field_functions_.size();
      ofile
        << "data" << ff << "=np.zeros([" << interpolation_points_.size()
        << ",5])\n";
    }

    offset = std::string("");
  }
  else if (Chi::mpi.process_count>1)
  {

    if (Chi::mpi.location_id != (Chi::mpi.process_count-1))
    {
      submod_name = base_name;
      submod_name = submod_name + std::to_string(Chi::mpi.location_id+1);

      ofile << "import " << submod_name << "\n\n";
    }


  }

  for (int ff=0; ff < field_functions_.size(); ff++)
  {
    const auto& ff_ctx = ff_contexts_[ff];

    if (Chi::mpi.process_count>1 and Chi::mpi.location_id!=0)
    {
      ofile
        << "def AddData" << ff << "(data" << ff << "):\n";

      offset = std::string("  ");
    }
    for (int p=0; p < interpolation_points_.size(); p++)
    {
      if ((not ff_ctx.interpolation_points_has_ass_cell[p])  &&
          (Chi::mpi.location_id != 0))
      {
        continue;
      }

      ofile << offset << "data" << ff << "[" << p << ",0] = "
            << interpolation_points_[p].x << "\n";
      ofile << offset << "data" << ff << "[" << p << ",1] = "
            << interpolation_points_[p].y << "\n";
      ofile << offset << "data" << ff << "[" << p << ",2] = "
            << interpolation_points_[p].z << "\n";

      double d = delta_d_ * p;

      ofile << offset << "data" << ff << "[" << p << ",3] = "
            << d << "\n";
      ofile << offset << "data" << ff << "[" << p << ",4] = "
            << ff_ctx.interpolation_points_values[p] << "\n";


    }

    ofile << offset << "done=True\n";
    ofile << "\n\n";
    if ((Chi::mpi.process_count>1) &&
        (Chi::mpi.location_id != (Chi::mpi.process_count-1)))
    {
      ofile << offset << submod_name
      << ".AddData" << ff << "(data" << ff << ")\n";
    }


  }

  for (int ca=0; ca < custom_arrays_.size(); ca++)
  {
    int ff = ca + field_functions_.size();

    if (Chi::mpi.process_count>1 and Chi::mpi.location_id!=0)
    {
      ofile
        << "def AddData" << ff << "(data" << ff << "):\n";

      offset = std::string("  ");
    }

    std::string op("= ");
    if (Chi::mpi.location_id != 0) op = std::string("+= ");

    for (int p=0; p < interpolation_points_.size(); p++)
    {
      ofile << offset << "data" << ff << "[" << p << ",0] = "
            << interpolation_points_[p].x << "\n";
      ofile << offset << "data" << ff << "[" << p << ",1] = "
            << interpolation_points_[p].y << "\n";
      ofile << offset << "data" << ff << "[" << p << ",2] = "
            << interpolation_points_[p].z << "\n";

      double d = delta_d_ * p;
      double value = 0.0;

      if (p < custom_arrays_[ca].size())
        value = custom_arrays_[ca][p];

      ofile << offset << "data" << ff << "[" << p << ",3] = "
            << d << "\n";
      ofile << offset << "data" << ff << "[" << p << ",4] " << op
            << value << "\n";
    }
    ofile << offset << "done=True\n";
    ofile << "\n\n";
    if ((Chi::mpi.process_count>1) &&
        (Chi::mpi.location_id != (Chi::mpi.process_count-1)))
    {
      ofile << offset << submod_name
            << ".AddData" << ff << "(data" << ff << ")\n";
    }
  }


  if (Chi::mpi.location_id == 0)
  {
    ofile << "plt.figure(1)\n";
    for (int ff=0; ff < field_functions_.size(); ff++)
    {
      ofile << "plt.plot(data" << ff << "[:,3],data"
            << ff << "[:,4]"
            << ",label=\"" << field_functions_[ff]->TextName() << "\""
            << ")\n";
    }
    for (int ca=0; ca < custom_arrays_.size(); ca++)
    {
      int ff = ca + field_functions_.size();
      ofile << "plt.plot(data" << ff << "[:,3],data"
            << ff << "[:,4]"
            << ",label=\"CustomArray" << ca << "\""
            << ")\n";
    }
    ofile << "plt.legend()\n"
             "plt.grid(which='major')\n";
    ofile << "plt.show()\n";
  }


  ofile.close();

  Chi::log.Log()
    << "Exported Python files for field func \""
    << field_functions_[0]->TextName()
    << "\" to base name \""
    << base_name << "\" Successfully";


}