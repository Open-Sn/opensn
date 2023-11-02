#ifndef CHI_PHYSICS_MATERIAL_PROPERTY_BASE_H
#define CHI_PHYSICS_MATERIAL_PROPERTY_BASE_H

#include <string>
#include <vector>

#include "chi_lua.h"

namespace chi_physics
{
enum class PropertyType
{
  SCALAR_VALUE = 1,
  TRANSPORT_XSECTIONS = 10,
  ISOTROPIC_MG_SOURCE = 11
};

//###################################################################
/** Base class for material properties.*/
class MaterialProperty
{
private:
  const PropertyType type_;

public:
  std::string property_name;

  explicit MaterialProperty(PropertyType in_type) : type_(in_type) {}

  virtual ~MaterialProperty() = default;

  PropertyType Type() { return type_; }

  virtual double GetScalarValue() { return 0.0; }

  virtual void PushLuaTable(lua_State* L) const;
};

} // namespace chi_physics

#endif // CHI_PHYSICS_MATERIAL_PROPERTY_BASE_H