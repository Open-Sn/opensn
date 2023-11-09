#pragma once

#include "framework/object.h"

namespace chi_physics::field_operations
{

/**The base field operation class.*/
class FieldOperation : public ChiObject
{
public:
  /**Returns the input parameters.*/
  static chi::InputParameters GetInputParameters();

  /**Constructor.*/
  explicit FieldOperation(const chi::InputParameters& params);

  virtual void Execute() = 0;

  virtual ~FieldOperation() = default;
};

} // namespace chi_physics::field_operations
