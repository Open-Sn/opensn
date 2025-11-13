// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/math/functions/function.h"

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
public:
  FieldFunctionInterpolationVolume()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::VOLUME),
      logical_volume_(nullptr),
      op_type_(FieldFunctionInterpolationOperation::OP_SUM),
      op_value_(0.0)
  {
  }

  ~FieldFunctionInterpolationVolume() override = default;

  std::shared_ptr<LogicalVolume> GetLogicalVolume() const { return logical_volume_; }

  void SetLogicalVolume(std::shared_ptr<LogicalVolume> lv) { logical_volume_ = lv; }

  FieldFunctionInterpolationOperation GetOperationType() const { return op_type_; }

  void SetOperationType(FieldFunctionInterpolationOperation op_type) { op_type_ = op_type; }

  double& GetOpValue() { return op_value_; }

  double GetValue() const { return op_value_; }

  void SetOperationFunction(const ScalarMaterialFunction& function) { oper_function_ = function; }

  void Initialize() override;

  void Execute() override;

private:
  std::shared_ptr<LogicalVolume> logical_volume_;
  FieldFunctionInterpolationOperation op_type_;
  double op_value_;
  ScalarMaterialFunction oper_function_;
  std::vector<std::uint32_t> cell_local_ids_inside_logvol_;

public:
  static std::shared_ptr<FieldFunctionInterpolationVolume> Create();
};

} // namespace opensn
