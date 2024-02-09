#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
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

OpenSnRegisterObjectInNamespace(physics, FieldFunctionGridBased);

InputParameters
FieldFunctionGridBased::GetInputParameters()
{
  InputParameters params = FieldFunction::GetInputParameters();

  params.SetDocGroup("DocFieldFunction");

  params.AddRequiredParameter<std::string>("sdm_type",
                                           "The spatial discretization type to be used");

  params.AddOptionalParameter(
    "initial_value", 0.0, "The initial value to assign to the field function");

  params.AddOptionalParameter("quadrature_order",
                              0,
                              "If supplied, will overwrite the default for the "
                              "specific discretization-coordinate system combination.");

  params.AddOptionalParameter(
    "coordinate_system", "cartesian", "Coordinate system to apply to element mappings");

  params.ConstrainParameterRange("sdm_type", AllowableRangeList::New({"FV", "PWLC", "PWLD"}));
  params.ConstrainParameterRange(
    "coordinate_system", AllowableRangeList::New({"cartesian", "cylindrical", "spherical"}));

  return params;
}

FieldFunctionGridBased::FieldFunctionGridBased(const InputParameters& params)
  : FieldFunction(params),
    sdm_(MakeSpatialDiscretization(params)),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(GetCurrentHandler().GetGrid()->GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(params.GetParamValue<double>("initial_value"));
}

FieldFunctionGridBased::FieldFunctionGridBased(
  const std::string& text_name,
  std::shared_ptr<SpatialDiscretization>& discretization_ptr,
  opensn::Unknown unknown)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(discretization_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
}

FieldFunctionGridBased::FieldFunctionGridBased(const std::string& text_name,
                                               std::shared_ptr<SpatialDiscretization>& sdm_ptr,
                                               opensn::Unknown unknown,
                                               const std::vector<double>& field_vector)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  ChiInvalidArgumentIf(field_vector.size() != ghosted_field_vector_->LocalSize(),
                       "Constructor called with incompatible size field vector.");

  ghosted_field_vector_->Set(field_vector);
}

FieldFunctionGridBased::FieldFunctionGridBased(const std::string& text_name,
                                               std::shared_ptr<SpatialDiscretization>& sdm_ptr,
                                               opensn::Unknown unknown,
                                               double field_value)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(field_value);
}

const SpatialDiscretization&
FieldFunctionGridBased::GetSpatialDiscretization() const
{
  return *sdm_;
}

const std::vector<double>&
FieldFunctionGridBased::FieldVectorRead() const
{
  return ghosted_field_vector_->LocalSTLData();
}

std::vector<double>&
FieldFunctionGridBased::FieldVector()
{
  return ghosted_field_vector_->LocalSTLData();
}

std::shared_ptr<SpatialDiscretization>
FieldFunctionGridBased::MakeSpatialDiscretization(const InputParameters& params)
{
  const auto& user_params = params.ParametersAtAssignment();
  const auto& grid_ptr = GetCurrentHandler().GetGrid();
  const auto sdm_type = params.GetParamValue<std::string>("sdm_type");

  typedef FiniteVolume FV;
  typedef PieceWiseLinearContinuous PWLC;
  typedef PieceWiseLinearDiscontinuous PWLD;

  if (sdm_type == "FV") return FV::New(*grid_ptr);

  CoordinateSystemType cs_type = CoordinateSystemType::CARTESIAN;
  std::string cs = "cartesian";
  if (user_params.Has("coordinate_system"))
  {
    cs = params.GetParamValue<std::string>("coordinate_system");

    if (cs == "cartesian") cs_type = CoordinateSystemType::CARTESIAN;
    if (cs == "cylindrical") cs_type = CoordinateSystemType::CYLINDRICAL;
    if (cs == "spherical") cs_type = CoordinateSystemType::SPHERICAL;
  }

  QuadratureOrder q_order = QuadratureOrder::SECOND;

  if (user_params.Has("quadrature_order"))
  {
    const uint32_t max_order = static_cast<uint32_t>(QuadratureOrder::FORTYTHIRD);
    const uint32_t q_order_int = params.GetParamValue<uint32_t>("quadrature_order");
    ChiInvalidArgumentIf(q_order_int > max_order,
                         "Invalid quadrature point order " + std::to_string(q_order_int));
    q_order = static_cast<QuadratureOrder>(q_order_int);
  }
  else // Defaulted
  {
    if (cs == "cartesian") q_order = QuadratureOrder::SECOND;
    if (cs == "cylindrical") q_order = QuadratureOrder::THIRD;
    if (cs == "spherical") q_order = QuadratureOrder::FOURTH;
  }

  // clang-format off
  if (sdm_type == "PWLC")
    return PWLC::New(*grid_ptr, q_order, cs_type);
  else if (sdm_type == "PWLD")
    return PWLD::New(*grid_ptr, q_order, cs_type);
  // clang-format on

  // If not returned by now
  ChiInvalidArgument("Unsupported sdm_type \"" + sdm_type + "\"");
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

void
FieldFunctionGridBased::UpdateFieldVector(const std::vector<double>& field_vector)
{
  ChiInvalidArgumentIf(field_vector.size() < ghosted_field_vector_->LocalSize(),
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

void
FieldFunctionGridBased::ExportMultipleToVTK(
  const std::string& file_base_name,
  const std::vector<std::shared_ptr<const FieldFunctionGridBased>>& ff_list)
{
  const std::string fname = "chi_physics::FieldFunction::ExportMultipleToVTK";
  log.Log() << "Exporting field functions to VTK with file base \"" << file_base_name << "\"";

  if (ff_list.empty())
    throw std::logic_error(fname + ": Cannot be used with empty field-function"
                                   " list");

  // Setup master and check slaves
  const auto& master_ff_ptr = ff_list.front();
  const auto& master_ff = *master_ff_ptr;

  for (const auto& ff_ptr : ff_list)
    if (ff_ptr != master_ff_ptr)
      if (&ff_ptr->sdm_->Grid() != &master_ff_ptr->sdm_->Grid())
        throw std::logic_error(fname +
                               ": Cannot be used with field functions based on different grids.");

  // Get grid
  const auto& grid = master_ff.sdm_->Grid();

  auto ugrid = PrepareVtkUnstructuredGrid(grid);

  // Upload cell/point data
  auto cell_data = ugrid->GetCellData();
  auto point_data = ugrid->GetPointData();
  for (const auto& ff_ptr : ff_list)
  {
    const auto field_vector = ff_ptr->GetGhostedFieldVector();

    const auto& uk_man = ff_ptr->GetUnknownManager();
    const auto& unknown = ff_ptr->Unknown();
    const auto& sdm = ff_ptr->sdm_;
    const size_t num_comps = unknown.NumComponents();

    for (uint c = 0; c < num_comps; ++c)
    {
      std::string component_name = ff_ptr->TextName() + unknown.text_name_;
      if (num_comps > 1) component_name += unknown.component_text_names_[c];

      vtkNew<vtkDoubleArray> point_array;
      vtkNew<vtkDoubleArray> cell_array;

      point_array->SetName(component_name.c_str());
      cell_array->SetName(component_name.c_str());

      // Populate the array here
      for (const auto& cell : grid.local_cells)
      {
        const size_t num_nodes = sdm->GetCellNumNodes(cell);

        if (num_nodes == cell.vertex_ids_.size())
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
          for (int n = 0; n < cell.vertex_ids_.size(); ++n)
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
}

std::vector<double>
FieldFunctionGridBased::GetGhostedFieldVector() const
{
  return ghosted_field_vector_->LocalSTLData();
}

std::vector<double>
FieldFunctionGridBased::GetPointValue(const Vector3& point) const
{
  typedef const int64_t cint64_t;
  const auto& uk_man = GetUnknownManager();
  const size_t num_components = uk_man.GetTotalUnknownStructureSize();

  size_t local_num_point_hits = 0;
  std::vector<double> local_point_value(num_components, 0.0);

  const auto& xyz_min = local_grid_bounding_box_.first;
  const auto& xyz_max = local_grid_bounding_box_.second;

  const double xmin = xyz_min.x;
  const double ymin = xyz_min.y;
  const double zmin = xyz_min.z;

  const double xmax = xyz_max.x;
  const double ymax = xyz_max.y;
  const double zmax = xyz_max.z;

  const auto& field_vector = *ghosted_field_vector_;

  if (point.x >= xmin and point.x <= xmax and point.y >= ymin and point.y <= ymax and
      point.z >= zmin and point.z <= zmax)
  {
    const auto& grid = sdm_->Grid();
    for (const auto& cell : grid.local_cells)
    {
      if (grid.CheckPointInsideCell(cell, point))
      {
        const auto& cell_mapping = sdm_->GetCellMapping(cell);
        std::vector<double> shape_values;
        cell_mapping.ShapeValues(point, shape_values);

        local_num_point_hits += 1;

        const size_t num_nodes = cell_mapping.NumNodes();
        for (size_t c = 0; c < num_components; ++c)
        {
          for (size_t j = 0; j < num_nodes; ++j)
          {
            cint64_t dof_map_j = sdm_->MapDOFLocal(cell, j, uk_man, 0, c);
            const double dof_value_j = field_vector[dof_map_j];

            local_point_value[c] += dof_value_j * shape_values[j];
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
FieldFunctionGridBased::Evaluate(const Cell& cell,
                                 const Vector3& position,
                                 unsigned int component) const
{
  const auto& field_vector = *ghosted_field_vector_;

  typedef const int64_t cint64_t;
  const auto& cell_mapping = sdm_->GetCellMapping(cell);

  std::vector<double> shape_values;
  cell_mapping.ShapeValues(position, shape_values);

  double value = 0.0;
  const size_t num_nodes = cell_mapping.NumNodes();
  for (size_t j = 0; j < num_nodes; ++j)
  {
    cint64_t dof_map = sdm_->MapDOFLocal(cell, j, GetUnknownManager(), 0, component);

    value += field_vector[dof_map] * shape_values[j];
  }

  return value;
}

} // namespace opensn
