#include "framework/materials/material_property.h"

#include "framework/object_factory.h"

namespace opensn
{

RegisterChiObject(chi_objects, MaterialProperty);

InputParameters
MaterialProperty::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>("name", "Text name associated with this property");

  return params;
}

MaterialProperty::MaterialProperty(const InputParameters& params)
  : name_(params.GetParamValue<std::string>("name"))
{
}

const std::string&
MaterialProperty::TextName() const
{
  return name_;
}

} // namespace opensn
