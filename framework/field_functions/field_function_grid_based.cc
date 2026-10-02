// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/mesh/mesh/grid_vtk_utils.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include "framework/runtime.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <petsc.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

namespace opensn
{

namespace
{

/**
 * Uniform bins over the bounding box of the local cells, for locating points in cells.
 *
 * Each bin lists, in increasing order, the local cells whose bounding box overlaps it. A cell
 * contains only points inside its bounding box, so the cells that can contain a point are the
 * ones listed in the bin of the point. The bounding boxes are padded so that this also holds for
 * points that Mesh::CheckPointInsideCell accepts within its tolerance. The bins are sized for
 * about one cell each on average, so building them and locating the points take time that grows
 * with the number of cells plus the number of points, instead of their product, unless the cell
 * sizes vary so strongly that some bins hold many cells.
 */
class LocalCellBins
{
public:
  LocalCellBins(const Mesh& grid, const Vector3& box_min, const Vector3& box_max) : origin_(box_min)
  {
    const auto& cells = grid.GetLocalCells();
    const std::array<double, 3> extent = {
      box_max.x - box_min.x, box_max.y - box_min.y, box_max.z - box_min.z};
    // The pad also exceeds the tolerance of Mesh::CheckPointInsideCell, which grows with the
    // magnitude of the coordinates.
    const double scale = std::max({extent[0], extent[1], extent[2], 1.0});
    const double magnitude = std::max({1.0,
                                       std::abs(box_min.x),
                                       std::abs(box_min.y),
                                       std::abs(box_min.z),
                                       std::abs(box_max.x),
                                       std::abs(box_max.y),
                                       std::abs(box_max.z)});
    pad_ = std::max(1.0e-9 * scale, 128.0 * std::numeric_limits<double>::epsilon() * magnitude);

    // Bins per axis, about one cell per bin, on the axes along which the mesh has extent. An axis
    // shorter than one bin gets a single bin and is dropped, and the bins are resized over the
    // other axes. The product of the bin counts of the remaining axes is then at most the number
    // of cells, so a thin domain does not get many more bins than cells.
    const auto num_cells = static_cast<double>(std::max<size_t>(cells.size(), 1));
    std::array<bool, 3> binned = {};
    for (int a = 0; a < 3; ++a)
      binned[a] = extent[a] > pad_;
    double bins_per_length = 0.0;
    for (bool dropped = true; dropped;)
    {
      double binned_volume = 1.0;
      int num_binned = 0;
      for (int a = 0; a < 3; ++a)
        if (binned[a])
        {
          binned_volume *= extent[a];
          ++num_binned;
        }
      bins_per_length =
        num_binned > 0 ? std::pow(num_cells / binned_volume, 1.0 / num_binned) : 0.0;
      dropped = false;
      for (int a = 0; a < 3; ++a)
        if (binned[a] and extent[a] * bins_per_length < 1.0)
        {
          binned[a] = false;
          dropped = true;
        }
    }
    for (int a = 0; a < 3; ++a)
    {
      num_bins_[a] = 1;
      if (binned[a])
        num_bins_[a] =
          static_cast<size_t>(std::clamp(std::floor(extent[a] * bins_per_length), 1.0, num_cells));
      width_[a] = extent[a] > pad_ ? extent[a] / static_cast<double>(num_bins_[a]) : 1.0;
    }

    // Two passes over the cells: count the cells of each bin, then list them.
    std::vector<std::array<size_t, 6>> ranges(cells.size());
    offsets_.assign(num_bins_[0] * num_bins_[1] * num_bins_[2] + 1, 0);
    for (size_t c = 0; c < cells.size(); ++c)
    {
      Vector3 lo = grid.GlobalVertex(cells[c].vertex_ids.front());
      Vector3 hi = lo;
      for (const auto vid : cells[c].vertex_ids)
      {
        const auto& v = grid.GlobalVertex(vid);
        lo = Vector3(std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z));
        hi = Vector3(std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z));
      }
      const std::array<double, 3> lo_a = {lo.x, lo.y, lo.z};
      const std::array<double, 3> hi_a = {hi.x, hi.y, hi.z};
      auto& range = ranges[c];
      for (int a = 0; a < 3; ++a)
      {
        range[a] = BinIndex(a, lo_a[a] - pad_);
        range[3 + a] = BinIndex(a, hi_a[a] + pad_);
      }
      ForEachBin(range, [this](size_t bin) { ++offsets_[bin + 1]; });
    }
    for (size_t b = 1; b < offsets_.size(); ++b)
      offsets_[b] += offsets_[b - 1];
    cell_indices_.resize(offsets_.back());
    std::vector<size_t> fill(offsets_.begin(), offsets_.end() - 1);
    for (size_t c = 0; c < cells.size(); ++c)
      ForEachBin(ranges[c], [&](size_t bin) { cell_indices_[fill[bin]++] = c; });
  }

  /// Return the indices of the local cells that can contain the point, in increasing order.
  std::span<const size_t> Candidates(const Vector3& point) const
  {
    if (not(std::isfinite(point.x) and std::isfinite(point.y) and std::isfinite(point.z)))
      return {};
    const size_t bin = (BinIndex(2, point.z) * num_bins_[1] + BinIndex(1, point.y)) * num_bins_[0] +
                       BinIndex(0, point.x);
    return {cell_indices_.data() + offsets_[bin], offsets_[bin + 1] - offsets_[bin]};
  }

private:
  size_t BinIndex(int axis, double coordinate) const
  {
    const std::array<double, 3> origin = {origin_.x, origin_.y, origin_.z};
    const double index = std::floor((coordinate - origin[axis]) / width_[axis]);
    return static_cast<size_t>(std::clamp(index, 0.0, static_cast<double>(num_bins_[axis] - 1)));
  }

  template <typename Function>
  void ForEachBin(const std::array<size_t, 6>& range, const Function& function) const
  {
    for (size_t k = range[2]; k <= range[5]; ++k)
      for (size_t j = range[1]; j <= range[4]; ++j)
        for (size_t i = range[0]; i <= range[3]; ++i)
          function((k * num_bins_[1] + j) * num_bins_[0] + i);
  }

  Vector3 origin_;
  double pad_ = 0.0;
  std::array<size_t, 3> num_bins_ = {1, 1, 1};
  std::array<double, 3> width_ = {1.0, 1.0, 1.0};
  std::vector<size_t> offsets_;
  std::vector<size_t> cell_indices_;
};

} // namespace

InputParameters
FieldFunctionGridBased::GetInputParameters()
{
  InputParameters params = FieldFunction::GetInputParameters();

  params.AddOptionalParameter(
    "discretization", "PWLC", "The spatial discretization type to be used");
  params.AddOptionalParameter(
    "coord_sys", "cartesian", "Coordinate system to apply to element mappings");
  params.AddOptionalParameter("quadrature_order",
                              0,
                              "If supplied, will overwrite the default for the "
                              "specific discretization-coordinate system combination.");

  params.AddOptionalParameter(
    "initial_value", 0.0, "The initial value to assign to the field function");

  params.ConstrainParameterRange("discretization", AllowableRangeList::New({"PWLC", "PWLD"}));
  params.ConstrainParameterRange(
    "coord_sys", AllowableRangeList::New({"cartesian", "cylindrical", "spherical"}));

  return params;
}

FieldFunctionGridBased::FieldFunctionGridBased(const InputParameters& params)
  : FieldFunction(params),
    discretization_(MakeSpatialDiscretization(params)),
    ghosted_field_vector_(MakeFieldVector(*discretization_, GetUnknownManager())),
    local_grid_bounding_box_(params.GetSharedPtrParam<Mesh>("mesh")->GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(params.GetParamValue<double>("initial_value"));
}

FieldFunctionGridBased::FieldFunctionGridBased(
  std::string name, std::shared_ptr<SpatialDiscretization>& discretization_ptr, Unknown unknown)
  : FieldFunction(std::move(name), std::move(unknown)),
    discretization_(discretization_ptr),
    ghosted_field_vector_(MakeFieldVector(*discretization_, GetUnknownManager())),
    local_grid_bounding_box_(discretization_->GetMesh()->GetLocalBoundingBox())
{
}

FieldFunctionGridBased::FieldFunctionGridBased(
  std::string name,
  std::shared_ptr<SpatialDiscretization>& discretization_ptr,
  Unknown unknown,
  const std::vector<double>& field_vector)
  : FieldFunction(std::move(name), std::move(unknown)),
    discretization_(discretization_ptr),
    ghosted_field_vector_(MakeFieldVector(*discretization_, GetUnknownManager())),
    local_grid_bounding_box_(discretization_->GetMesh()->GetLocalBoundingBox())
{
  OpenSnInvalidArgumentIf(field_vector.size() != ghosted_field_vector_->GetLocalSize(),
                          "Constructor called with incompatible size field vector.");

  ghosted_field_vector_->Set(field_vector);
}

FieldFunctionGridBased::FieldFunctionGridBased(
  std::string name,
  std::shared_ptr<SpatialDiscretization>& discretization_ptr,
  Unknown unknown,
  double field_value)
  : FieldFunction(std::move(name), std::move(unknown)),
    discretization_(discretization_ptr),
    ghosted_field_vector_(MakeFieldVector(*discretization_, GetUnknownManager())),
    local_grid_bounding_box_(discretization_->GetMesh()->GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(field_value);
}

const SpatialDiscretization&
FieldFunctionGridBased::GetSpatialDiscretization() const
{
  return *discretization_;
}

std::vector<double>&
FieldFunctionGridBased::GetLocalFieldVector()
{
  return ghosted_field_vector_->GetLocalSTLData();
}

const std::vector<double>&
FieldFunctionGridBased::GetLocalFieldVector() const
{
  return ghosted_field_vector_->GetLocalSTLData();
}

std::vector<double>
FieldFunctionGridBased::GetGhostedFieldVector() const
{
  return ghosted_field_vector_->GetLocalSTLData();
}

void
FieldFunctionGridBased::UpdateFieldVector(const std::vector<double>& field_vector)
{
  OpenSnInvalidArgumentIf(field_vector.size() < ghosted_field_vector_->GetLocalSize(),
                          "Attempted update with a vector of insufficient size.");

  ghosted_field_vector_->Set(field_vector);
  ghosted_field_vector_->CommunicateGhostEntries();
}

void
FieldFunctionGridBased::UpdateFieldVector(const Vec& field_vector)
{
  ghosted_field_vector_->CopyLocalValues(field_vector);

  ghosted_field_vector_->CommunicateGhostEntries();
}

bool
FieldFunctionGridBased::CanUpdate() const
{
  return static_cast<bool>(update_callback_) and
         (not can_update_callback_ or can_update_callback_());
}

void
FieldFunctionGridBased::Update()
{
  OpenSnLogicalErrorIf(not CanUpdate(),
                       "Field function \"" + GetName() + "\" is not currently updateable.");
  update_callback_(*this);
}

void
FieldFunctionGridBased::SetUpdateCallback(UpdateCallback callback,
                                          CanUpdateCallback can_update_callback)
{
  update_callback_ = std::move(callback);
  can_update_callback_ = std::move(can_update_callback);
}

std::vector<double>
FieldFunctionGridBased::GetPointValue(const Vector3& point) const
{
  const auto& uk_man = GetUnknownManager();
  const size_t num_components = uk_man.GetTotalUnknownStructureSize();

  size_t local_num_point_hits = 0;
  std::vector<double> local_point_value(num_components, 0.0);

  const auto& xyz_min = local_grid_bounding_box_.first;
  const auto& xyz_max = local_grid_bounding_box_.second;

  const auto xmin = xyz_min.x;
  const auto ymin = xyz_min.y;
  const auto zmin = xyz_min.z;

  const auto xmax = xyz_max.x;
  const auto ymax = xyz_max.y;
  const auto zmax = xyz_max.z;

  const auto& field_vector = *ghosted_field_vector_;

  if (point.x >= xmin and point.x <= xmax and point.y >= ymin and point.y <= ymax and
      point.z >= zmin and point.z <= zmax)
  {
    const auto& grid = discretization_->GetMesh();
    for (std::uint32_t cell_local_id = 0; cell_local_id < grid->GetLocalCellCount();
         ++cell_local_id)
    {
      const auto& cell = grid->GetLocalCell(cell_local_id);
      if (grid->CheckPointInsideCell(cell, point))
      {
        const auto& cell_mapping = discretization_->GetLocalCellMapping(cell_local_id);
        Vector<double> shape_values;
        cell_mapping.ShapeValues(point, shape_values);

        local_num_point_hits += 1;

        const auto num_nodes = cell_mapping.GetNumNodes();
        for (size_t c = 0; c < num_components; ++c)
        {
          for (size_t j = 0; j < num_nodes; ++j)
          {
            const auto dof_map = discretization_->MapDOFLocal(cell_local_id, j, uk_man, 0, c);
            const double dof_value = field_vector[dof_map];

            local_point_value[c] += dof_value * shape_values(j);
          } // for node i
        } // for component c
      } // if inside cell
    } // for cell
  } // if in bounding box

  // Communicate number of point hits
  size_t global_num_point_hits = 0;
  mpi_comm.all_reduce(local_num_point_hits, global_num_point_hits, mpi::op::sum<size_t>());
  OpenSnLogicalErrorIf(global_num_point_hits == 0,
                       "FieldFunctionGridBased::GetPointValue: No cell identified containing "
                       "the specified point.");

  std::vector<double> global_point_value(num_components, 0.0);
  OpenSnLogicalErrorIf(num_components > static_cast<size_t>(std::numeric_limits<int>::max()),
                       "FieldFunctionGridBased::GetPointValue: Number of components exceeds "
                       "MPI reduction count limit.");
  const auto point_value_size = static_cast<int>(num_components);
  mpi_comm.all_reduce(
    local_point_value.data(), point_value_size, global_point_value.data(), mpi::op::sum<double>());

  Scale(global_point_value, 1.0 / static_cast<double>(global_num_point_hits));

  return global_point_value;
}

std::vector<std::vector<std::vector<double>>>
FieldFunctionGridBased::GetPointValues(const FFList& field_functions,
                                       const std::vector<Vector3>& points)
{
  const std::string fname = "FieldFunctionGridBased::GetPointValues";
  OpenSnInvalidArgumentIf(field_functions.empty(), fname + ": Field-function list is empty.");
  OpenSnInvalidArgumentIf(not field_functions.front(),
                          fname + ": Field-function list contains null.");

  const auto& master = *field_functions.front();
  const auto master_grid = master.discretization_->GetMesh();
  for (const auto& field_function : field_functions)
  {
    OpenSnInvalidArgumentIf(not field_function, fname + ": Field-function list contains null.");
    OpenSnInvalidArgumentIf(field_function->discretization_->GetMesh() != master_grid,
                            fname + ": Field functions use different meshes.");
  }

  std::vector<size_t> component_offsets(field_functions.size() + 1, 0);
  for (size_t f = 0; f < field_functions.size(); ++f)
    component_offsets[f + 1] =
      component_offsets[f] + field_functions[f]->GetUnknownManager().GetTotalUnknownStructureSize();
  const size_t num_components = component_offsets.back();

  OpenSnInvalidArgumentIf(num_components != 0 and
                            points.size() > std::numeric_limits<size_t>::max() / num_components,
                          fname + ": Result size exceeds addressable memory.");
  const size_t num_values = points.size() * num_components;
  OpenSnInvalidArgumentIf(points.size() > static_cast<size_t>(std::numeric_limits<int>::max()) or
                            num_values > static_cast<size_t>(std::numeric_limits<int>::max()),
                          fname + ": MPI reduction count limit exceeded.");

  std::vector<size_t> local_num_point_hits(points.size(), 0);
  std::vector<double> local_values(num_values, 0.0);
  const auto& [xyz_min, xyz_max] = master.local_grid_bounding_box_;
  const auto& local_cells = master_grid->GetLocalCells();
  const LocalCellBins bins(*master_grid, xyz_min, xyz_max);

  for (size_t p = 0; p < points.size(); ++p)
  {
    const auto& point = points[p];
    if (point.x < xyz_min.x or point.x > xyz_max.x or point.y < xyz_min.y or point.y > xyz_max.y or
        point.z < xyz_min.z or point.z > xyz_max.z)
      continue;

    for (const size_t cell_local_id : bins.Candidates(point))
    {
      const auto& cell = local_cells[cell_local_id];
      if (not master_grid->CheckPointInsideCell(cell, point))
        continue;

      ++local_num_point_hits[p];
      const auto& master_mapping = master.discretization_->GetLocalCellMapping(cell_local_id);
      Vector<double> master_shape_values;
      master_mapping.ShapeValues(point, master_shape_values);
      for (size_t f = 0; f < field_functions.size(); ++f)
      {
        const auto& field_function = *field_functions[f];
        const auto& discretization = *field_function.discretization_;
        const auto& unknown_manager = field_function.GetUnknownManager();
        const auto& cell_mapping = discretization.GetLocalCellMapping(cell_local_id);
        const auto& field_vector = *field_function.ghosted_field_vector_;
        Vector<double> alternate_shape_values;
        const Vector<double>* shape_values = &master_shape_values;
        if (field_function.discretization_ != master.discretization_)
        {
          cell_mapping.ShapeValues(point, alternate_shape_values);
          shape_values = &alternate_shape_values;
        }

        const size_t num_nodes = cell_mapping.GetNumNodes();
        const size_t field_num_components = unknown_manager.GetTotalUnknownStructureSize();
        for (size_t c = 0; c < field_num_components; ++c)
          for (size_t j = 0; j < num_nodes; ++j)
          {
            const auto dof_map =
              discretization.MapDOFLocal(cell_local_id, j, unknown_manager, 0, c);
            const size_t value_index = p * num_components + component_offsets[f] + c;
            local_values[value_index] += field_vector[dof_map] * (*shape_values)(j);
          }
      }
    }
  }

  std::vector<size_t> global_num_point_hits(points.size(), 0);
  if (not points.empty())
    mpi_comm.all_reduce(local_num_point_hits.data(),
                        static_cast<int>(points.size()),
                        global_num_point_hits.data(),
                        mpi::op::sum<size_t>());
  for (size_t p = 0; p < points.size(); ++p)
    OpenSnLogicalErrorIf(global_num_point_hits[p] == 0,
                         fname + ": No cell identified containing point " + std::to_string(p) +
                           ".");

  std::vector<double> global_values(num_values, 0.0);
  if (not global_values.empty())
    mpi_comm.all_reduce(local_values.data(),
                        static_cast<int>(num_values),
                        global_values.data(),
                        mpi::op::sum<double>());

  std::vector<std::vector<std::vector<double>>> result(field_functions.size());
  for (size_t f = 0; f < field_functions.size(); ++f)
  {
    const size_t field_num_components = component_offsets[f + 1] - component_offsets[f];
    result[f].assign(points.size(), std::vector<double>(field_num_components, 0.0));
    for (size_t p = 0; p < points.size(); ++p)
      for (size_t c = 0; c < field_num_components; ++c)
        result[f][p][c] = global_values[p * num_components + component_offsets[f] + c] /
                          static_cast<double>(global_num_point_hits[p]);
  }
  return result;
}

double
FieldFunctionGridBased::Evaluate(const Cell& cell, const Vector3& position, int component) const
{
  const auto& field_vector = *ghosted_field_vector_;
  const auto mesh = discretization_->GetMesh();
  const auto cell_local_id = mesh->MapCellGlobalID2LocalID(cell.global_id);
  const auto& cell_mapping = discretization_->GetLocalCellMapping(cell_local_id);

  Vector<double> shape_values;
  cell_mapping.ShapeValues(position, shape_values);

  double value = 0.0;
  const size_t num_nodes = cell_mapping.GetNumNodes();
  for (size_t j = 0; j < num_nodes; ++j)
  {
    const auto dof_map =
      discretization_->MapDOFLocal(cell_local_id, j, GetUnknownManager(), 0, component);
    value += field_vector[dof_map] * shape_values(j);
  }

  return value;
}

void
FieldFunctionGridBased::ExportMultipleToPVTU(
  const std::string& file_base_name,
  const std::vector<std::shared_ptr<const FieldFunctionGridBased>>& ff_list)
{
  const std::string fname = "FieldFunctionGridBased::ExportMultipleToPVTU";
  log.Log() << "Exporting field functions to PVTU with file base \"" << file_base_name << "\"";

  if (ff_list.empty())
    throw std::logic_error(fname + ": Cannot be used with empty field-function"
                                   " list");

  // Setup master and check slaves
  const auto& master_ff_ptr = ff_list.front();
  const auto& master_ff = *master_ff_ptr;

  for (const auto& ff_ptr : ff_list)
    if (ff_ptr != master_ff_ptr)
      if (ff_ptr->discretization_->GetMesh() != master_ff_ptr->discretization_->GetMesh())
        throw std::logic_error(fname +
                               ": Cannot be used with field functions based on different grids.");

  // Get grid
  const auto& grid = master_ff.discretization_->GetMesh();

  auto ugrid = PrepareVtkUnstructuredGrid(grid);

  // Upload cell/point data
  auto* point_data = ugrid->GetPointData();
  for (const auto& ff_ptr : ff_list)
  {
    const auto field_vector = ff_ptr->GetGhostedFieldVector();

    const auto& uk_man = ff_ptr->GetUnknownManager();
    const auto& unknown = ff_ptr->GetUnknown();
    const auto& sdm = ff_ptr->discretization_;
    const size_t num_comps = unknown.GetNumComponents();

    for (uint c = 0; c < num_comps; ++c)
    {
      std::string component_name = ff_ptr->GetName() + unknown.name;
      if (num_comps > 1)
        component_name += unknown.component_names[c];

      vtkNew<vtkDoubleArray> point_array;

      point_array->SetName(component_name.c_str());

      // Populate the array here
      for (std::uint32_t cell_local_id = 0; cell_local_id < grid->GetLocalCellCount();
           ++cell_local_id)
      {
        const auto& cell = grid->GetLocalCell(cell_local_id);
        const size_t num_nodes = sdm->GetCellNumNodes(cell);

        if (num_nodes == cell.vertex_ids.size())
        {
          for (size_t n = 0; n < num_nodes; ++n)
          {
            const auto nmap = sdm->MapDOFLocal(cell_local_id, n, uk_man, 0, c);

            const double field_value = field_vector[nmap];

            point_array->InsertNextValue(field_value);
          } // for node
        }
        else
        {
          double node_average = 0.0;
          for (size_t n = 0; n < num_nodes; ++n)
          {
            const auto nmap = sdm->MapDOFLocal(cell_local_id, n, uk_man, 0, c);

            const double field_value = field_vector[nmap];
            node_average += field_value;
          } // for node
          node_average /= static_cast<double>(num_nodes);
          for (std::size_t n = 0; n < cell.vertex_ids.size(); ++n)
          {
            point_array->InsertNextValue(node_average);
          } // for vertex
        }

      } // for cell

      point_data->AddArray(point_array);
    } // for component
  } // for ff_ptr

  WritePVTUFiles(ugrid, file_base_name);

  log.Log() << "Done exporting field functions to PVTU.";
  opensn::mpi_comm.barrier();
}

std::shared_ptr<SpatialDiscretization>
FieldFunctionGridBased::MakeSpatialDiscretization(const InputParameters& params)
{
  const auto grid = params.GetSharedPtrParam<Mesh>("mesh");
  const auto sdm_type = params.GetParamValue<std::string>("discretization");

  std::string cs = "cartesian";
  if (params.IsParameterValid("coord_sys"))
    cs = params.GetParamValue<std::string>("coord_sys");

  QuadratureOrder q_order = QuadratureOrder::SECOND;

  if (params.IsParameterValid("quadrature_order"))
  {
    const auto max_order = static_cast<uint32_t>(QuadratureOrder::FORTYTHIRD);
    const auto q_order_int = params.GetParamValue<uint32_t>("quadrature_order");
    OpenSnInvalidArgumentIf(q_order_int > max_order,
                            "Invalid quadrature point order " + std::to_string(q_order_int));
    q_order = static_cast<QuadratureOrder>(q_order_int);
  }
  else // Defaulted
  {
    if (cs == "cartesian")
      q_order = QuadratureOrder::SECOND;
    else if (cs == "cylindrical")
      q_order = QuadratureOrder::THIRD;
    else if (cs == "spherical")
      q_order = QuadratureOrder::FOURTH;
    else
      throw std::logic_error("Unrecognized coordinate system.");
  }

  if (sdm_type == "PWLC")
    return PieceWiseLinearContinuous::New(grid, q_order);
  else if (sdm_type == "PWLD")
    return PieceWiseLinearDiscontinuous::New(grid, q_order);

  // If not returned by now
  OpenSnInvalidArgument("Unsupported discretization \"" + sdm_type + "\"");
}

std::unique_ptr<GhostedParallelSTLVector>
FieldFunctionGridBased::MakeFieldVector(const SpatialDiscretization& discretization,
                                        const UnknownManager& uk_man)
{
  auto field = std::make_unique<GhostedParallelSTLVector>(discretization.GetNumLocalDOFs(uk_man),
                                                          discretization.GetNumGlobalDOFs(uk_man),
                                                          discretization.GetGhostDOFIndices(uk_man),
                                                          mpi_comm);

  return field;
}

} // namespace opensn
