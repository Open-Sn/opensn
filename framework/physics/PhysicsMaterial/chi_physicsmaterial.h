#pragma once

#include "opensn/framework/physics/chi_physics_namespace.h"
#include "opensn/framework/physics/PhysicsMaterial/material_property_base.h"

#include <vector>
#include <memory>

namespace chi_physics
{

//###################################################################
/** Base class for materials used in physics simulations.*/
class Material
{
public:
  std::vector<std::shared_ptr<MaterialProperty>> properties_{};
  std::string name_ = "Unnamed Material";
};

} // namespace chi_physics
