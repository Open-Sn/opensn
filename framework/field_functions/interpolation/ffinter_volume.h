// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/math/functions/scalar_material_function.h"

namespace opensn
{

/**
 * Volume-wise field function interpolation.
 *
 * This interpolator allows the user to obtain quantities by logical
 * volume. If no logical volume is assigned to the method it will
 * default to operating over the entire volume.
 */
class FieldFunctionInterpolationVolume : public FieldFunctionInterpolation
{
private:
  std::shared_ptr<LogicalVolume> logical_volume_;
  FieldFunctionInterpolationOperation op_type_;
  double op_value_;
  std::shared_ptr<ScalarMaterialFunction> oper_function_;
  std::vector<uint64_t> cell_local_ids_inside_logvol_;

public:
  FieldFunctionInterpolationVolume()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::VOLUME),
      logical_volume_(nullptr),
      op_type_(FieldFunctionInterpolationOperation::OP_SUM),
      op_value_(0.0)
  {}

  virtual ~FieldFunctionInterpolationVolume() {}

  std::shared_ptr<LogicalVolume>& GetLogicalVolume() { return logical_volume_; }

  FieldFunctionInterpolationOperation& GetOperationType() { return op_type_; }

  double& GetOpValue() { return op_value_; }

  void SetOperationFunction(std::shared_ptr<ScalarMaterialFunction> function);

  void Initialize() override;

  void Execute() override;
};

} // namespace opensn
