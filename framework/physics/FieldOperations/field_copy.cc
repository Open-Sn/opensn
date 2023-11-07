#include "framework/physics/FieldOperations/field_copy.h"

#include "framework/ChiObjectFactory.h"
#include "framework/math/SpatialDiscretization/SpatialDiscretization.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_physics::field_operations
{

RegisterChiObject(chi_physics::field_operations, FieldCopyOperation);

/**Returns the input parameters.*/
chi::InputParameters
FieldCopyOperation::GetInputParameters()
{
  chi::InputParameters params = FieldOperation::GetInputParameters();

  params.SetDocGroup("DocFieldOperation");

  params.AddRequiredParameter<size_t>("to",
                                      "Handle to the field function that should "
                                      "receive the result of the operation.");

  params.AddRequiredParameter<size_t>("from",
                                      "Handle to the field function that should "
                                      "be copied from.");

  params.AddOptionalParameterArray(
    "to_components",
    std::vector<size_t>{},
    "List of component numbers that need to "
    "be written to. If this parameter is supplied then"
    "\"from_components\" also need to be specified (and vice versa) and they "
    "need to be the same length. If neither are supplied, all components are "
    "copied.");

  params.AddOptionalParameterArray(
    "from_components",
    std::vector<size_t>{},
    "List of component numbers that need to "
    "read from. If this parameter is supplied then"
    "\"to_components\" also need to be specified (and vice versa) and they "
    "need to be the same length. If neither are supplied, all components are "
    "copied.");

  return params;
}

/**Constructor.*/
FieldCopyOperation::FieldCopyOperation(const chi::InputParameters& params)
  : FieldOperation(params),
    to_field_handle_(params.GetParamValue<size_t>("to")),
    from_field_handle_(params.GetParamValue<size_t>("from"))
{
  // Get field functions
  {
    auto to_base_ptr =
      Chi::GetStackItemPtr(Chi::field_function_stack, to_field_handle_, __FUNCTION__);

    to_ff_ = std::dynamic_pointer_cast<FieldFunctionGridBased>(to_base_ptr);

    ChiLogicalErrorIf(not to_ff_,
                      "\"to\" field function must be based on "
                      "FieldFunctionGridBased");
  }

  {
    auto from_base_ptr =
      Chi::GetStackItemPtr(Chi::field_function_stack, from_field_handle_, __FUNCTION__);

    from_ff_ = std::dynamic_pointer_cast<FieldFunctionGridBased>(from_base_ptr);

    ChiLogicalErrorIf(not to_ff_,
                      "\"to\" field function must be based on "
                      "FieldFunctionGridBased");
  }

  // Check number of components are compatible
  const auto& user_supplied_params = params.ParametersAtAssignment();

  if (user_supplied_params.Has("to_components") and
      (not user_supplied_params.Has("from_components")))
    ChiInvalidArgument("If \"to_components\" is specified then "
                       "\"from_components\" must also be specified");

  if (user_supplied_params.Has("from_components") and
      (not user_supplied_params.Has("to_components")))
    ChiInvalidArgument("If \"from_components\" is specified then "
                       "\"to_components\" must also be specified");

  if (user_supplied_params.Has("to_components") and user_supplied_params.Has("from_components"))
  {
    to_components_ = user_supplied_params.GetParamVectorValue<size_t>("to_components");
    from_components_ = user_supplied_params.GetParamVectorValue<size_t>("from_components");

    ChiInvalidArgumentIf(to_components_.size() != from_components_.size(),
                         "\"to_components\" and \"from_components\" must have "
                         "the same number of entries");
  }
  else
  {
    ChiInvalidArgumentIf(to_ff_->GetUnknownManager().GetTotalUnknownStructureSize() !=
                           from_ff_->GetUnknownManager().GetTotalUnknownStructureSize(),
                         "The number of components of the unknowns in the field functions are"
                         " not compatible");

    const size_t num_comps = to_ff_->GetUnknownManager().GetTotalUnknownStructureSize();
    to_components_.reserve(num_comps);
    from_components_.reserve(num_comps);
    for (size_t c = 0; c < num_comps; ++c)
    {
      to_components_.push_back(c);
      from_components_.push_back(c);
    }
  }

  // Check grids are compatible
  ChiInvalidArgumentIf(std::addressof(to_ff_->GetSpatialDiscretization().Grid()) !=
                         std::addressof(from_ff_->GetSpatialDiscretization().Grid()),
                       "Currently the two field functions must operate on the "
                       "same grid");
}

void
FieldCopyOperation::Execute()
{
  typedef const int64_t cint64_t;
  const auto& sdm = to_ff_->GetSpatialDiscretization();
  const auto& uk_man = to_ff_->GetUnknownManager();
  const auto& grid = sdm.Grid();

  const size_t num_comps = to_components_.size();

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto nodes_xyz = cell_mapping.GetNodeLocations();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t c = 0; c < num_comps; ++c)
      {
        const size_t cto = to_components_[c];
        const size_t cfrom = from_components_[c];

        const double value = from_ff_->Evaluate(cell, nodes_xyz[i], cfrom);

        cint64_t dof_map = sdm.MapDOFLocal(cell, i, uk_man, 0, cto);

        to_ff_->FieldVector()[dof_map] = value;
      } // for component c
    }   // for node i
  }     // for cell
}

} // namespace chi_physics::field_operations
