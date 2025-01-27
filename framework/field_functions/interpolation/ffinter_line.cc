// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/interpolation/ffinter_line.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/vector_ghost_communicator/vector_ghost_communicator.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/mesh/cell/cell.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <fstream>

namespace opensn
{

std::shared_ptr<FieldFunctionInterpolationLine>
FieldFunctionInterpolationLine::Create()
{
  auto ffi = std::make_shared<FieldFunctionInterpolationLine>();
  field_func_interpolation_stack.emplace_back(ffi);
  return ffi;
}

void
FieldFunctionInterpolationLine::Initialize()
{
  log.Log0Verbose1() << "Initializing line interpolator.";

  // Check for empty FF-list
  if (field_functions_.empty())
    throw std::logic_error("Unassigned field function in line field function interpolator.");

  // Create points;
  const Vector3 vif = pf_ - pi_;
  double delta_d = vif.Norm() / (number_of_points_ - 1);
  const auto omega = vif.Normalized();
  std::vector<Vector3> tmp_points(number_of_points_);
  tmp_points[0] = pi_;
  for (int k = 1; k < number_of_points_; ++k)
    tmp_points[k] = pi_ + omega * delta_d * k;

  ref_ff_ = field_functions_.front();
  const auto& sdm = ref_ff_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();

  // Find local points and associated cells
  auto estimated_local_size = number_of_points_ / opensn::mpi_comm.size();
  local_interpolation_points_.reserve(estimated_local_size);
  local_cells_.reserve(estimated_local_size);
  for (const auto& cell : grid.local_cells)
  {
    for (int p = 0; p < number_of_points_; ++p)
    {
      auto& point = tmp_points[p];
      if (grid.CheckPointInsideCell(cell, point))
      {
        local_interpolation_points_.push_back(point);
        local_cells_.push_back(cell.local_id);
      }
    }
  }

  log.Log0Verbose1() << "Finished initializing interpolator.";
}

void
FieldFunctionInterpolationLine::Execute()
{
  log.Log0Verbose1() << "Executing line interpolator.";

  const auto& sdm = ref_ff_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& uk_man = ref_ff_->GetUnknownManager();
  const auto uid = 0;
  const auto cid = ref_component_;
  const auto field_data = ref_ff_->GetGhostedFieldVector();

  double local_max = 0.0, local_sum = 0.0, local_avg = 0.0;
  size_t local_size = local_interpolation_points_.size();
  local_interpolation_values_.resize(local_size);
  for (auto p = 0; p < local_size; ++p)
  {
    auto& point = local_interpolation_points_[p];
    auto cell_local_index = local_cells_[p];
    const auto& cell = grid.local_cells[cell_local_index];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    Vector<double> shape_function_vals(num_nodes, 0.0);
    cell_mapping.ShapeValues(point, shape_function_vals);
    double point_value = 0.0;
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);
      point_value += shape_function_vals(i) * field_data[imap];
    }
    local_interpolation_values_[p] = point_value;
    local_max = std::max(point_value, local_max);
    local_sum += point_value;
  }

  if (op_type_ == FieldFunctionInterpolationOperation::OP_SUM)
    mpi_comm.all_reduce(local_sum, op_value_, mpi::op::sum<double>());
  else if (op_type_ == FieldFunctionInterpolationOperation::OP_AVG)
  {
    size_t global_size = 0;
    mpi_comm.all_reduce(local_size, global_size, mpi::op::sum<size_t>());
    double global_sum = 0.0;
    mpi_comm.all_reduce(local_sum, global_sum, mpi::op::sum<double>());
    op_value_ = global_sum / static_cast<double>(global_size);
  }
  else if (op_type_ == FieldFunctionInterpolationOperation::OP_MAX)
    mpi_comm.all_reduce(local_max, op_value_, mpi::op::max<double>());
}

void
FieldFunctionInterpolationLine::ExportToCSV(std::string base_name) const
{
  // Populate local coordinate and interpolation data
  std::vector<double> local_data;
  local_data.reserve(4 * number_of_points_ / opensn::mpi_comm.size());
  for (auto p = 0; p < local_interpolation_points_.size(); ++p)
  {
    auto& point = local_interpolation_points_[p];
    local_data.push_back(point.x);
    local_data.push_back(point.y);
    local_data.push_back(point.z);
    local_data.push_back(local_interpolation_values_[p]);
  }

  // Compute size of local coordinate and interpolation data and send to rank 0
  std::vector<int> local_data_sizes(opensn::mpi_comm.size(), 0);
  int local_size = local_data.size();
  opensn::mpi_comm.gather(local_size, local_data_sizes, 0);

  // Compute global size of coordinate and interpolation data and each location's offset
  std::vector<double> global_data;
  std::vector<int> offsets(opensn::mpi_comm.size(), 0);
  int global_data_size = 0;
  if (opensn::mpi_comm.rank() == 0)
  {
    offsets.resize(opensn::mpi_comm.size());
    for (int i = 0; i < opensn::mpi_comm.size(); ++i)
    {
      offsets[i] = global_data_size;
      global_data_size += local_data_sizes[i];
    }
    global_data.resize(global_data_size);
  }

  // Send all local data to rank 0
  opensn::mpi_comm.gather(local_data, global_data, local_data_sizes, offsets, 0);

  if (opensn::mpi_comm.rank() == 0)
  {
    // Zip data and sort data by coordinate
    std::vector<std::tuple<double, double, double, double>> values;
    values.reserve(global_data_size / 4);
    for (size_t i = 0; i < global_data_size; i += 4)
    {
      values.emplace_back(std::make_tuple(
        global_data[i], global_data[i + 1], global_data[i + 2], global_data[i + 3]));
    }

    std::stable_sort(values.begin(),
                     values.end(),
                     [](const std::tuple<double, double, double, double>& a,
                        const std::tuple<double, double, double, double>& b)
                     { return std::get<2>(a) < std::get<2>(b); });
    std::stable_sort(values.begin(),
                     values.end(),
                     [](const std::tuple<double, double, double, double>& a,
                        const std::tuple<double, double, double, double>& b)
                     { return std::get<1>(a) < std::get<1>(b); });
    std::stable_sort(values.begin(),
                     values.end(),
                     [](const std::tuple<double, double, double, double>& a,
                        const std::tuple<double, double, double, double>& b)
                     { return std::get<0>(a) < std::get<0>(b); });

    // Write sorted data to CSV file
    std::ofstream ofile;
    std::string filename = base_name + "_" + ref_ff_->GetName() + std::string(".csv");
    ofile.open(filename);
    ofile << "x,y,z," << ref_ff_->GetName() << "\n";
    for (auto point_data : values)
    {
      auto [x, y, z, fv] = point_data;
      ofile << std::setprecision(7) << x << "," << y << "," << z << "," << fv << "\n";
    }
    ofile.close();

    log.Log() << "Exported CSV file for field func \"" << ref_ff_->GetName() << "\" to \""
              << filename << "\"";
  }
}

} // namespace opensn
