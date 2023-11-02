#ifndef CHITECH_CHI_MATERIAL_PROPERTY_H
#define CHITECH_CHI_MATERIAL_PROPERTY_H

#include "ChiObject.h"

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

}

#endif // CHITECH_CHI_MATERIAL_PROPERTY_H
