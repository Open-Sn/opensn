#ifndef CHITECH_FIELD_OPERATION_H
#define CHITECH_FIELD_OPERATION_H

#include "ChiObject.h"

namespace chi_physics::field_operations
{

/**The base field operation class.*/
class FieldOperation : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();

  explicit FieldOperation(const chi::InputParameters& params);

  virtual void Execute() = 0;

  virtual ~FieldOperation() = default;
};

} // namespace chi_physics::field_operations

#endif // CHITECH_FIELD_OPERATION_H
