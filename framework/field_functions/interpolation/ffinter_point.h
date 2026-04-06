// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/data_types/vector3.h"

namespace opensn
{

/// A point based interpolation function.
class FieldFunctionInterpolationPoint : public FieldFunctionInterpolation
{
public:
  FieldFunctionInterpolationPoint()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::POINT),
      locally_owned_(false),
      owning_cell_gid_(0),
      point_value_(0.0)
  {
  }

  ~FieldFunctionInterpolationPoint() override = default;

  Vector3 GetPointOfInterest() const { return point_of_interest_; }

  void SetPointOfInterest(const Vector3& point) { point_of_interest_ = point; }

  void Execute() override;

  /// Gets the value of the field function evaluation at the point.
  double GetPointValue() const;

protected:
  void RebuildPointLocationData();

  Vector3 point_of_interest_;
  bool locally_owned_;
  uint64_t owning_cell_gid_;
  double point_value_;

public:
  static std::shared_ptr<FieldFunctionInterpolationPoint> Create();
};

} // namespace opensn
