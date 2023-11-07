#pragma once

#include "framework/physics/field_operations/field_operation.h"
#include "framework/physics/field_function/field_function_grid_based.h"

namespace chi_physics::field_operations
{

/**Field operaiton that copies components of one field to the
 * components of another.*/
class FieldCopyOperation : public FieldOperation
{
private:
  const size_t to_field_handle_;
  const size_t from_field_handle_;

  std::vector<size_t> to_components_;
  std::vector<size_t> from_components_;

  std::shared_ptr<FieldFunctionGridBased> to_ff_;
  std::shared_ptr<const FieldFunctionGridBased> from_ff_;

public:
  /**Returns the input parameters.*/
  static chi::InputParameters GetInputParameters();

  /**Constructor.*/
  explicit FieldCopyOperation(const chi::InputParameters& params);

  void Execute() override;
};

} // namespace chi_physics::field_operations
