#include "framework/physics/FieldFunction/fieldfunction.h"

#include "framework/logging/chi_log_exceptions.h"

namespace chi_physics
{

// ##################################################################
/**Returns required input parameters.*/
chi::InputParameters
FieldFunction::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>("name",
                                           "Named to be associated with this field function");

  params.AddOptionalParameter(
    "unknown_type", "Scalar", "The type of the variable for this field function");

  params.AddOptionalParameter("num_components",
                              1,
                              "The number of components to attach to the variable. "
                              "Only effective when \"type\" is VectorN.");

  //=============== Constrain values
  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "unknown_type", AllowableRangeList::New({"Scalar", "Vector2", "Vector3", "VectorN"}));

  params.ConstrainParameterRange("num_components", AllowableRangeLowLimit::New(1));

  return params;
}

// ##################################################################
/**ObjectMaker based constructor.*/
FieldFunction::FieldFunction(const chi::InputParameters& params)
  : ChiObject(params),
    text_name_(params.GetParamValue<std::string>("name")),
    unknown_((params.GetParamValue<std::string>("unknown_type") == "Scalar")
               ? chi_math::Unknown(chi_math::UnknownType::SCALAR)
             : (params.GetParamValue<std::string>("unknown_type") == "Vector2")
               ? chi_math::Unknown(chi_math::UnknownType::VECTOR_2)
             : (params.GetParamValue<std::string>("unknown_type") == "Vector3")
               ? chi_math::Unknown(chi_math::UnknownType::VECTOR_2)
             : (params.GetParamValue<std::string>("unknown_type") == "VectorN")
               ? chi_math::Unknown(chi_math::UnknownType::VECTOR_N,
                                   params.GetParamValue<unsigned int>("num_components"))
               : chi_math::Unknown(chi_math::UnknownType::SCALAR)),
    unknown_manager_({unknown_})
{
}

// ##################################################################
/**Conventional constructor.*/
FieldFunction::FieldFunction(const std::string& text_name, chi_math::Unknown unknown)
  : text_name_(text_name), unknown_(std::move(unknown)), unknown_manager_({unknown_})
{
}

// ##################################################################
/**Stack change to `chi::field_function_stack.*/
void
FieldFunction::PushOntoStack(std::shared_ptr<ChiObject>& new_object)
{
  auto ff_ptr = std::dynamic_pointer_cast<FieldFunction>(new_object);

  ChiLogicalErrorIf(not ff_ptr, "Bad trouble when casting object to field function");

  Chi::field_function_stack.push_back(ff_ptr);
  new_object->SetStackID(Chi::field_function_stack.size() - 1);
}

} // namespace chi_physics
