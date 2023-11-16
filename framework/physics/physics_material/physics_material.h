#pragma once

#include "framework/physics/physics_namespace.h"
#include "framework/physics/physics_material/material_property_base.h"

#include <vector>
#include <memory>

namespace opensn
{

/** Base class for materials used in physics simulations.*/
class Material
{
public:
  std::vector<std::shared_ptr<MaterialProperty>> properties_{};
  std::string name_ = "Unnamed Material";
};

} // namespace opensn
