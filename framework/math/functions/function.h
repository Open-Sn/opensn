#pragma once

#include "framework/object.h"

namespace opensn
{

/**
 * Base class for functions
 */
class Function : public Object
{
public:
  static InputParameters GetInputParameters();

protected:
  explicit Function(const InputParameters& params);
};

} // namespace opensn
