// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/mesh/mesh.h"

namespace opensn
{

/** A line based interpolation function.*/
class FieldFunctionInterpolationLine : public FieldFunctionInterpolation
{
private:
  int number_of_points_;
  double op_value_;
  Vector3 pi_, pf_;
  std::vector<uint64_t> local_cells_;
  std::vector<Vector3> local_interpolation_points_;
  std::vector<double> local_interpolation_values_;
  std::shared_ptr<FieldFunctionGridBased> ref_ff_;
  FieldFunctionInterpolationOperation op_type_;

public:
  FieldFunctionInterpolationLine()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::LINE),
      number_of_points_(2),
      op_value_(0.0),
      op_type_(FieldFunctionInterpolationOperation::OP_SUM)
  {
  }

  virtual ~FieldFunctionInterpolationLine() {}

  int& GetNumberOfPoints() { return number_of_points_; }

  void SetNumberOfPoints(int number) { number_of_points_ = number; }

  Vector3 GetInitialPoint() { return pi_; }

  void SetInitialPoint(Vector3 point) { pi_ = point; }

  Vector3 GetFinalPoint() { return pf_; }

  void SetFinalPoint(Vector3 point) { pf_ = point; }

  void Initialize() override;

  void Execute() override;

  FieldFunctionInterpolationOperation& GetOperationType() { return op_type_; }

  void SetOperationType(FieldFunctionInterpolationOperation type) { op_type_ = type; }

  double GetOpValue() { return op_value_; }

  std::string GetDefaultFileBaseName() const override { return "ZLFFI"; }

  void ExportPython(std::string base_name) override;
};

} // namespace opensn
