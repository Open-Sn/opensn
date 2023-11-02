#ifndef CHITECH_FIELDFUNCTIONINTERFACE_H
#define CHITECH_FIELDFUNCTIONINTERFACE_H

#include "parameters/input_parameters.h"

namespace chi_physics
{

class FieldFunction;

/**Interface class to add a dependency on a logical volume. Two things need to
* be done to use this interface. 1) Derive from it. 2) Add its parameters to
* the child class. Now it will require a handle to a FieldFunction in
* the input language.*/
class FieldFunctionInterface
{
protected:
  static chi::InputParameters GetInputParameters();

  explicit FieldFunctionInterface(const chi::InputParameters& params);

  FieldFunction* GetFieldFunction() const;

private:
  chi::ParameterBlock field_function_param_;
};

}

#endif // CHITECH_FIELDFUNCTIONINTERFACE_H
