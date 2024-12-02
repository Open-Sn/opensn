// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_continuum/grid_vtk_utils.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <petsc.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

namespace opensn
{

InputParameters
FieldFunctionGridBased::GetInputParameters()
{
  InputParameters params = FieldFunction::GetInputParameters();

  params.SetDocGroup("DocFieldFunction");

  params.AddOptionalParameter("discretization", "FV", "The spatial discretization type to be used");
  params.AddOptionalParameter(
    "coordinate_system", "cartesian", "Coordinate system to apply to element mappings");
  params.AddOptionalParameter("quadrature_order",
                              0,
                              "If supplied, will overwrite the default for the "
                              "specific discretization-coordinate system combination.");

  params.AddOptionalParameter(
    "initial_value", 0.0, "The initial value to assign to the field function");

  params.ConstrainParameterRange("discretization", AllowableRangeList::New({"FV", "PWLC", "PWLD"}));
  params.ConstrainParameterRange(
    "coordinate_system", AllowableRangeList::New({"cartesian", "cylindrical", "spherical"}));

  return params;
}

FieldFunctionGridBased::FieldFunctionGridBased(const InputParameters& params)
  : FieldFunction(params),
    discretization_(MakeSpatialDiscretization(params)),
    ghosted_field_vector_(MakeFieldVector(*discretization_, GetUnknownManager())),
    local_grid_bounding_box_(GetCurrentMesh()->GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(params.GetParamValue<double>("initial_value"));
}

FieldFunctionGridBased::FieldFunctionGridBased(
  std::string name, std::shared_ptr<SpatialDiscretization>& discretization_ptr, Unknown unknown)
  : FieldFunction(std::move(name), std::move(unknown)),
    discretization_(discretization_ptr),
    ghosted_field_vector_(MakeFieldVector(*discretization_, GetUnknownManager())),
    local_grid_bounding_box_(discretization_->Grid().GetLocalBoundingBox())
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
    local_grid_bounding_box_(discretization_->Grid().GetLocalBoundingBox())
{
  OpenSnInvalidArgumentIf(field_vector.size() != ghosted_field_vector_->LocalSize(),
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
    local_grid_bounding_box_(discretization_->Grid().GetLocalBoundingBox())
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
  return ghosted_field_vector_->LocalSTLData();
}

const std::vector<double>&
FieldFunctionGridBased::GetLocalFieldVector() const
{
  return ghosted_field_vector_->LocalSTLData();
}

std::vector<double>
FieldFunctionGridBased::GetGhostedFieldVector() const
{
  return ghosted_field_vector_->LocalSTLData();
}

void
FieldFunctionGridBased::UpdateFieldVector(const std::vector<double>& field_vector)
{
  OpenSnInvalidArgumentIf(field_vector.size() < ghosted_field_vector_->LocalSize(),
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
    const auto& grid = discretization_->Grid();
    for (const auto& cell : grid.local_cells)
    {
      if (grid.CheckPointInsideCell(cell, point))
      {
        const auto& cell_mapping = discretization_->GetCellMapping(cell);
        Vector<double> shape_values;
        cell_mapping.ShapeValues(point, shape_values);

        local_num_point_hits += 1;

        const auto num_nodes = cell_mapping.NumNodes();
        for (size_t c = 0; c < num_components; ++c)
        {
          for (size_t j = 0; j < num_nodes; ++j)
          {
            const auto dof_map = discretization_->MapDOFLocal(cell, j, uk_man, 0, c);
            const double dof_value = field_vector[dof_map];

            local_point_value[c] += dof_value * shape_values(j);
          } // for node i
        }   // for component c
      }     // if inside cell
    }       // for cell
  }         // if in bounding box

  // Communicate number of point hits
  size_t globl_num_point_hits;
  mpi_comm.all_reduce(local_num_point_hits, globl_num_point_hits, mpi::op::sum<size_t>());

  std::vector<double> globl_point_value(num_components, 0.0);
  mpi_comm.all_reduce(
    local_point_value.data(), 1, globl_point_value.data(), mpi::op::sum<double>());

  Scale(globl_point_value, 1.0 / static_cast<double>(globl_num_point_hits));

  return globl_point_value;
}

double
FieldFunctionGridBased::Evaluate(const Cell& cell, const Vector3& position, int component) const
{
  const auto& field_vector = *ghosted_field_vector_;

  const auto& cell_mapping = discretization_->GetCellMapping(cell);

  Vector<double> shape_values;
  cell_mapping.ShapeValues(position, shape_values);

  double value = 0.0;
  const size_t num_nodes = cell_mapping.NumNodes();
  for (size_t j = 0; j < num_nodes; ++j)
  {
    const auto dof_map = discretization_->MapDOFLocal(cell, j, GetUnknownManager(), 0, component);
    value += field_vector[dof_map] * shape_values(j);
  }

  return value;
}

void
FieldFunctionGridBased::ExportMultipleToVTK(
  const std::string& file_base_name,
  const std::vector<std::shared_ptr<const FieldFunctionGridBased>>& ff_list)
{
  const std::string fname = "FieldFunctionGridBased::ExportMultipleToVTK";
  log.Log() << "Exporting field functions to VTK with file base \"" << file_base_name << "\"";

  if (ff_list.empty())
    throw std::logic_error(fname + ": Cannot be used with empty field-function"
                                   " list");

  // Setup master and check slaves
  const auto& master_ff_ptr = ff_list.front();
  const auto& master_ff = *master_ff_ptr;

  for (const auto& ff_ptr : ff_list)
    if (ff_ptr != master_ff_ptr)
      if (&ff_ptr->discretization_->Grid() != &master_ff_ptr->discretization_->Grid())
        throw std::logic_error(fname +
                               ": Cannot be used with field functions based on different grids.");

  // Get grid
  const auto& grid = master_ff.discretization_->Grid();

  auto ugrid = PrepareVtkUnstructuredGrid(grid);

  // Upload cell/point data
  auto cell_data = ugrid->GetCellData();
  auto point_data = ugrid->GetPointData();
  for (const auto& ff_ptr : ff_list)
  {
    const auto field_vector = ff_ptr->GetGhostedFieldVector();

    const auto& uk_man = ff_ptr->GetUnknownManager();
    const auto& unknown = ff_ptr->GetUnknown();
    const auto& sdm = ff_ptr->discretization_;
    const size_t num_comps = unknown.NumComponents();

    for (uint c = 0; c < num_comps; ++c)
    {
      std::string component_name = ff_ptr->Name() + unknown.name;
      if (num_comps > 1)
        component_name += unknown.component_names[c];

      vtkNew<vtkDoubleArray> point_array;
      vtkNew<vtkDoubleArray> cell_array;

      point_array->SetName(component_name.c_str());
      cell_array->SetName(component_name.c_str());

      // Populate the array here
      for (const auto& cell : grid.local_cells)
      {
        const size_t num_nodes = sdm->GetCellNumNodes(cell);

        if (num_nodes == cell.vertex_ids.size())
        {
          double node_average = 0.0;
          for (int n = 0; n < num_nodes; ++n)
          {
            const int64_t nmap = sdm->MapDOFLocal(cell, n, uk_man, 0, c);

            const double field_value = field_vector[nmap];

            point_array->InsertNextValue(field_value);
            node_average += field_value;
          } // for node
          node_average /= static_cast<double>(num_nodes);
          cell_array->InsertNextValue(node_average);
        }
        else
        {
          double node_average = 0.0;
          for (int n = 0; n < num_nodes; ++n)
          {
            const int64_t nmap = sdm->MapDOFLocal(cell, n, uk_man, 0, c);

            const double field_value = field_vector[nmap];
            node_average += field_value;
          } // for node
          node_average /= static_cast<double>(num_nodes);
          cell_array->InsertNextValue(node_average);
          for (int n = 0; n < cell.vertex_ids.size(); ++n)
          {
            point_array->InsertNextValue(node_average);
          } // for vertex
        }

      } // for cell

      point_data->AddArray(point_array);
      cell_data->AddArray(cell_array);
    } // for component
  }   // for ff_ptr

  WritePVTUFiles(ugrid, file_base_name);

  log.Log() << "Done exporting field functions to VTK.";
  opensn::mpi_comm.barrier();
}

std::shared_ptr<SpatialDiscretization>
FieldFunctionGridBased::MakeSpatialDiscretization(const InputParameters& params)
{
  const auto& grid_ptr = GetCurrentMesh();
  const auto sdm_type = params.GetParamValue<std::string>("discretization");

  if (sdm_type == "FV")
    return FiniteVolume::New(*grid_ptr);

  CoordinateSystemType cs_type = CoordinateSystemType::CARTESIAN;
  std::string cs = "cartesian";
  if (params.IsParameterValid("coordinate_system"))
  {
    cs = params.GetParamValue<std::string>("coordinate_system");

    if (cs == "cartesian")
      cs_type = CoordinateSystemType::CARTESIAN;
    if (cs == "cylindrical")
      cs_type = CoordinateSystemType::CYLINDRICAL;
    if (cs == "spherical")
      cs_type = CoordinateSystemType::SPHERICAL;
  }

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
    if (cs == "cylindrical")
      q_order = QuadratureOrder::THIRD;
    if (cs == "spherical")
      q_order = QuadratureOrder::FOURTH;
  }

  if (sdm_type == "PWLC")
    return PieceWiseLinearContinuous::New(*grid_ptr, q_order, cs_type);
  else if (sdm_type == "PWLD")
    return PieceWiseLinearDiscontinuous::New(*grid_ptr, q_order, cs_type);

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
