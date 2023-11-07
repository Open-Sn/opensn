#pragma once

#include "framework/object.h"

namespace chi
{

/**Base class for a material property.*/
class MaterialProperty : public ChiObject
{
private:
  const std::string name_;

public:
  static chi::InputParameters GetInputParameters();
  explicit MaterialProperty(const chi::InputParameters& params);

  const std::string& TextName() const;
};

} // namespace chi
