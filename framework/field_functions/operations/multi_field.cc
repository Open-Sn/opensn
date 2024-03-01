#include "framework/field_functions/operations/multi_field.h"
#include "framework/object_factory.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(physics::field_operations, MultiFieldOperation);

InputParameters
MultiFieldOperation::GetInputParameters()
{
  InputParameters params = FieldOperation::GetInputParameters();

  params.SetDocGroup("DocFieldOperation");

  params.AddRequiredParameter<size_t>("result_field_handle",
                                      "Handle to the field function that should "
                                      "receive the result of the operation.");

  params.AddOptionalParameterArray(
    "dependent_field_handles", std::vector<int>{}, "List of dependent field handles");

  params.AddOptionalParameterArray("dependent_component_references",
                                   std::vector<unsigned int>{},
                                   "The components of interest for each "
                                   "handle supplied");
  params.AddOptionalParameterArray("result_component_references",
                                   std::vector<unsigned int>{},
                                   "The components to write into the result "
                                   "field function");

  params.AddRequiredParameter<size_t>("function_handle",
                                      "Handle to a FunctionDimAToDimB derived"
                                      " object");

  return params;
}

MultiFieldOperation::MultiFieldOperation(const InputParameters& params)
  : FieldOperation(params),
    result_field_handle_(params.GetParamValue<size_t>("result_field_handle")),
    dependent_field_handles_(params.GetParamVectorValue<size_t>("dependent_field_handles")),
    function_handle_(params.GetParamValue<size_t>("function_handle"))
{
  // Make component references
  const auto& user_supplied_params = params.ParametersAtAssignment();
  if (user_supplied_params.Has("dependent_component_references"))
  {
    dependent_field_ref_component_ =
      user_supplied_params.GetParamVectorValue<unsigned int>("dependent_component_references");

    OpenSnInvalidArgumentIf(dependent_field_ref_component_.size() !=
                              dependent_field_handles_.size(),
                            "Not each dependent field handle has an associated "
                            "reference component");
  }
  else
    dependent_field_ref_component_.assign(dependent_field_handles_.size(), 0);

  if (user_supplied_params.Has("result_component_references"))
  {
    result_component_references_ =
      user_supplied_params.GetParamVectorValue<unsigned int>("result_component_references");
  }
  else
    result_component_references_ = {0};

  // Process handles
  auto ff_base_ptr = GetStackItemPtr(field_function_stack, result_field_handle_, __FUNCTION__);

  primary_ff_ = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff_base_ptr);

  OpenSnLogicalErrorIf(not primary_ff_,
                       "Primary field function must be based on "
                       "FieldFunctionGridBased");

  for (const size_t dep_handle : dependent_field_handles_)
  {
    auto dep_ff_base_ptr = GetStackItemPtr(field_function_stack, dep_handle, __FUNCTION__);

    auto dep_ff_ptr = std::dynamic_pointer_cast<FieldFunctionGridBased>(dep_ff_base_ptr);

    OpenSnLogicalErrorIf(not dep_ff_ptr,
                         "Dependent field function must be based on "
                         "FieldFunctionGridBased");

    dependent_ffs_.push_back(dep_ff_ptr);
  }

  auto function_base_obj = GetStackItemPtr(object_stack, function_handle_, __FUNCTION__);

  function_ptr_ = std::dynamic_pointer_cast<FunctionDimAToDimB>(function_base_obj);

  OpenSnLogicalErrorIf(not function_ptr_, "Casting failure of function");

  const size_t num_dependencies = dependent_ffs_.size() + 4;
  OpenSnInvalidArgumentIf(function_ptr_->InputDimension() != num_dependencies,
                          std::string("The function will be called with "
                                      "x,y,z,material_id and then each of the "
                                      "dependent field function values. ") +
                            "This means the function needs to be callable"
                            " with 4+" +
                            std::to_string(dependent_ffs_.size()) +
                            " values, "
                            "however, it only supports " +
                            std::to_string(function_ptr_->InputDimension()) + " values.");

  OpenSnInvalidArgumentIf(function_ptr_->OutputDimension() != result_component_references_.size(),
                          std::string("The function will return ") +
                            std::to_string(function_ptr_->OutputDimension()) +
                            " values however this operation only requires " +
                            std::to_string(result_component_references_.size()) + " value(s).");
}

void
MultiFieldOperation::Execute()
{
  typedef unsigned int uint;
  typedef const int64_t cint64_t;
  const auto& sdm = primary_ff_->GetSpatialDiscretization();
  const auto& uk_man = primary_ff_->GetUnknownManager();
  const auto& grid = sdm.Grid();

  const size_t num_deps = dependent_ffs_.size();

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto nodes_xyz = cell_mapping.GetNodeLocations();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      std::vector<double> input_params;
      const auto& node_i_xyz = nodes_xyz[i];
      for (size_t d = 0; d < 3; ++d)
        input_params.push_back(node_i_xyz[d]);
      input_params.push_back(static_cast<double>(cell.material_id_));

      for (size_t k = 0; k < num_deps; ++k)
      {
        const auto& dep_ff = dependent_ffs_[k];
        const double value = dep_ff->Evaluate(cell, node_i_xyz, dependent_field_ref_component_[k]);
        input_params.push_back(value);
      }

      std::vector<double> output_params = function_ptr_->Evaluate(input_params);

      OpenSnLogicalErrorIf(output_params.size() != result_component_references_.size(),
                           "Calling function number of output values not matching");

      size_t k = 0;
      for (uint c : result_component_references_)
      {
        cint64_t dof_map = sdm.MapDOFLocal(cell, i, uk_man, 0, c);

        primary_ff_->FieldVector()[dof_map] = output_params[k++];
      }

    } // for node i
  }   // for cell
}

} // namespace opensn
