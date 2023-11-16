#include "framework/materials/material.h"

#include "framework/object_factory.h"

namespace opensn
{

RegisterChiObject(chi_objects, Material);

InputParameters
Material::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name", "The text name that will be associated with this material.");

  params.AddRequiredParameterArray(
    "properties", "Expects an array object handles that represents the properties.");

  return params;
}

Material::Material(const InputParameters& params) : name_(params.GetParamValue<std::string>("name"))
{
}

} // namespace opensn
