#pragma once

#include "ChiObject.h"

namespace chi
{

/**A generic material object used to group together multiple properties.*/
class Material : public ChiObject
{
private:
  std::string name_;

public:
  static InputParameters GetInputParameters();
  explicit Material(const chi::InputParameters& params);
};

} // namespace chi


