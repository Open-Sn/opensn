#pragma once

#include "framework/field_functions/field_function_interface.h"

namespace opensn
{

class FieldFunctionGridBased;

/**Interface class to add a dependency on a logical volume. Two things need to
 * be done to use this interface. 1) Derive from it. 2) Add its parameters to
 * the child class. Now it will require a handle to a GridBasedFieldFunction in
 * the input language.*/
class GridBasedFieldFunctionInterface : public FieldFunctionInterface
{
public:
  static InputParameters GetInputParameters();

  explicit GridBasedFieldFunctionInterface(const InputParameters& params);

  FieldFunctionGridBased* GetGridBasedFieldFunction() const;
};

} // namespace opensn
