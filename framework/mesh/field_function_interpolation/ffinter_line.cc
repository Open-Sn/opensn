#include "framework/mesh/field_function_interpolation/ffinter_line.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/vector_ghost_communicator//vector_ghost_communicator.h"
#include "framework/physics/field_function/field_function_grid_based.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi.h"
#include <fstream>

namespace chi_mesh
{

void
chi_mesh::FieldFunctionInterpolationLine::Initialize()
{
  Chi::log.Log0Verbose1() << "Initializing line interpolator.";
  // Check for empty FF-list
  if (field_functions_.empty())
    throw std::logic_error("Unassigned field function in line "
                           "field function interpolator.");

  // Create points;
  const chi_mesh::Vector3 vif = pf_ - pi_;
  delta_d_ = vif.Norm() / (number_of_points_ - 1);

  const auto omega = vif.Normalized();

  interpolation_points_.push_back(pi_);
  for (int k = 1; k < (number_of_points_); k++)
    interpolation_points_.push_back(pi_ + omega * delta_d_ * k);

  // Loop over contexts
  const size_t num_ff = field_functions_.size();
  for (size_t ff = 0; ff < num_ff; ff++)
  {
    ff_contexts_.emplace_back();
    auto& ff_context = ff_contexts_.back();

    ff_context.ref_ff = field_functions_[ff];
    const auto& sdm = ff_context.ref_ff->GetSpatialDiscretization();
    const auto& grid = sdm.Grid();

    ff_context.interpolation_points_ass_cell.assign(number_of_points_, 0);
    ff_context.interpolation_points_has_ass_cell.assign(number_of_points_, false);

    // Find a home for each point
    for (const auto& cell : grid.local_cells)
    {
      for (int p = 0; p < number_of_points_; p++)
      {
        const auto& point = interpolation_points_[p];
        if (grid.CheckPointInsideCell(cell, point))
        {
          ff_context.interpolation_points_ass_cell[p] = cell.local_id_;
          ff_context.interpolation_points_has_ass_cell[p] = true;
        }
      } // for point p
    }   // for cell
  }     // for ff

  Chi::log.Log0Verbose1() << "Finished initializing interpolator.";
}

void
FieldFunctionInterpolationLine::Execute()
{
  Chi::log.Log0Verbose1() << "Executing line interpolator.";
  for (int ff = 0; ff < field_functions_.size(); ff++)
  {
    auto& ff_ctx = ff_contexts_[ff];
    const auto& ref_ff = *ff_ctx.ref_ff;
    const auto& sdm = ref_ff.GetSpatialDiscretization();
    const auto& grid = sdm.Grid();

    const auto& uk_man = ref_ff.GetUnknownManager();
    const auto uid = 0;
    const auto cid = ref_component_;

    const auto field_data = ref_ff.GetGhostedFieldVector();

    ff_ctx.interpolation_points_values.assign(number_of_points_, 0.0);
    for (int p = 0; p < number_of_points_; ++p)
    {
      if (not ff_ctx.interpolation_points_has_ass_cell[p]) continue;

      const auto cell_local_index = ff_ctx.interpolation_points_ass_cell[p];
      const auto& cell = grid.local_cells[cell_local_index];
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      std::vector<double> shape_function_vals(num_nodes, 0.0);
      cell_mapping.ShapeValues(interpolation_points_[p], shape_function_vals);

      double point_value = 0.0;
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);

        point_value += shape_function_vals[i] * field_data[imap];
      } // for node i
      ff_ctx.interpolation_points_values[p] = point_value;
    } // for p
  }   // for ff
}

void
FieldFunctionInterpolationLine::ExportPython(std::string base_name)
{
  std::ofstream ofile;

  std::string fileName = base_name;
  fileName = fileName + std::to_string(Chi::mpi.location_id);
  fileName = fileName + std::string(".py");
  ofile.open(fileName);

  ofile << "import numpy as np\n"
           "import matplotlib.pyplot as plt\n"
        << "\n";

  std::string offset;
  std::string submod_name;
  if (Chi::mpi.location_id == 0)
  {
    submod_name = base_name;
    submod_name = submod_name + std::to_string(Chi::mpi.location_id + 1);

    if (Chi::mpi.process_count > 1) { ofile << "import " << submod_name << "\n\n"; }

    for (int ff = 0; ff < field_functions_.size(); ff++)
    {
      ofile << "data" << ff << "=np.zeros([" << interpolation_points_.size() << ",5])\n";
    }
    for (int ca = 0; ca < custom_arrays_.size(); ca++)
    {
      int ff = ca + field_functions_.size();
      ofile << "data" << ff << "=np.zeros([" << interpolation_points_.size() << ",5])\n";
    }

    offset = std::string("");
  }
  else if (Chi::mpi.process_count > 1)
  {

    if (Chi::mpi.location_id != (Chi::mpi.process_count - 1))
    {
      submod_name = base_name;
      submod_name = submod_name + std::to_string(Chi::mpi.location_id + 1);

      ofile << "import " << submod_name << "\n\n";
    }
  }

  for (int ff = 0; ff < field_functions_.size(); ff++)
  {
    const auto& ff_ctx = ff_contexts_[ff];

    if (Chi::mpi.process_count > 1 and Chi::mpi.location_id != 0)
    {
      ofile << "def AddData" << ff << "(data" << ff << "):\n";

      offset = std::string("  ");
    }
    for (int p = 0; p < interpolation_points_.size(); p++)
    {
      if ((not ff_ctx.interpolation_points_has_ass_cell[p]) && (Chi::mpi.location_id != 0))
      {
        continue;
      }

      ofile << offset << "data" << ff << "[" << p << ",0] = " << interpolation_points_[p].x << "\n";
      ofile << offset << "data" << ff << "[" << p << ",1] = " << interpolation_points_[p].y << "\n";
      ofile << offset << "data" << ff << "[" << p << ",2] = " << interpolation_points_[p].z << "\n";

      double d = delta_d_ * p;

      ofile << offset << "data" << ff << "[" << p << ",3] = " << d << "\n";
      ofile << offset << "data" << ff << "[" << p
            << ",4] = " << ff_ctx.interpolation_points_values[p] << "\n";
    }

    ofile << offset << "done=True\n";
    ofile << "\n\n";
    if ((Chi::mpi.process_count > 1) && (Chi::mpi.location_id != (Chi::mpi.process_count - 1)))
    {
      ofile << offset << submod_name << ".AddData" << ff << "(data" << ff << ")\n";
    }
  }

  for (int ca = 0; ca < custom_arrays_.size(); ca++)
  {
    int ff = ca + field_functions_.size();

    if (Chi::mpi.process_count > 1 and Chi::mpi.location_id != 0)
    {
      ofile << "def AddData" << ff << "(data" << ff << "):\n";

      offset = std::string("  ");
    }

    std::string op("= ");
    if (Chi::mpi.location_id != 0) op = std::string("+= ");

    for (int p = 0; p < interpolation_points_.size(); p++)
    {
      ofile << offset << "data" << ff << "[" << p << ",0] = " << interpolation_points_[p].x << "\n";
      ofile << offset << "data" << ff << "[" << p << ",1] = " << interpolation_points_[p].y << "\n";
      ofile << offset << "data" << ff << "[" << p << ",2] = " << interpolation_points_[p].z << "\n";

      double d = delta_d_ * p;
      double value = 0.0;

      if (p < custom_arrays_[ca].size()) value = custom_arrays_[ca][p];

      ofile << offset << "data" << ff << "[" << p << ",3] = " << d << "\n";
      ofile << offset << "data" << ff << "[" << p << ",4] " << op << value << "\n";
    }
    ofile << offset << "done=True\n";
    ofile << "\n\n";
    if ((Chi::mpi.process_count > 1) && (Chi::mpi.location_id != (Chi::mpi.process_count - 1)))
    {
      ofile << offset << submod_name << ".AddData" << ff << "(data" << ff << ")\n";
    }
  }

  if (Chi::mpi.location_id == 0)
  {
    ofile << "plt.figure(1)\n";
    for (int ff = 0; ff < field_functions_.size(); ff++)
    {
      ofile << "plt.plot(data" << ff << "[:,3],data" << ff << "[:,4]"
            << ",label=\"" << field_functions_[ff]->TextName() << "\""
            << ")\n";
    }
    for (int ca = 0; ca < custom_arrays_.size(); ca++)
    {
      int ff = ca + field_functions_.size();
      ofile << "plt.plot(data" << ff << "[:,3],data" << ff << "[:,4]"
            << ",label=\"CustomArray" << ca << "\""
            << ")\n";
    }
    ofile << "plt.legend()\n"
             "plt.grid(which='major')\n";
    ofile << "plt.show()\n";
  }

  ofile.close();

  Chi::log.Log() << "Exported Python files for field func \"" << field_functions_[0]->TextName()
                 << "\" to base name \"" << base_name << "\" Successfully";
}

} // namespace chi_mesh
