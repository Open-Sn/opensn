// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/operations/field_operation.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/functions/function_dimA_to_dimB.h"

namespace opensn
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

  std::shared_ptr<const FunctionDimAToDimB> function_ptr_;

public:
  /**Returns the input parameters.*/
  static InputParameters GetInputParameters();

  /**Constructor.*/
  explicit MultiFieldOperation(const InputParameters& params);

  void Execute() override;
};

} // namespace opensn
