#pragma once

#include "field_operation.h"
#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/Functions/function_dimA_to_dimB.h"

namespace chi_physics::field_operations
{

/**A field operation to manipulate a single field on the hand
 * of a number of other fields.*/
class MultiFieldOperation : public FieldOperation
{
private:
  const size_t result_field_handle_;
  const std::vector<size_t> dependent_field_handles_;
  const size_t function_handle_;

  std::vector<unsigned int> dependent_field_ref_component_;
  std::vector<unsigned int> result_component_references_;

  std::shared_ptr<FieldFunctionGridBased> primary_ff_;
  std::vector<std::shared_ptr<const FieldFunction>> dependent_ffs_;

  std::shared_ptr<const chi_math::FunctionDimAToDimB> function_ptr_;

public:
  static chi::InputParameters GetInputParameters();

  explicit MultiFieldOperation(const chi::InputParameters& params);

  void Execute() override;
};

} // namespace chi_physics::field_operations
