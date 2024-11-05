// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/mesh/mesh.h"

namespace opensn
{

/// A line based interpolation function.
class FieldFunctionInterpolationPoint : public FieldFunctionInterpolation
{
protected:
  Vector3 point_of_interest_;
  bool locally_owned_;
  uint64_t owning_cell_gid_;
  double point_value_;

public:
  FieldFunctionInterpolationPoint()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::POINT),
      locally_owned_(false),
      owning_cell_gid_(0),
      point_value_(0.0)
  {
  }

  virtual ~FieldFunctionInterpolationPoint() {}

  Vector3& PointOfInterest() { return point_of_interest_; }

  void Initialize() override;

  void Execute() override;

  /// Gets the value of the field function evaluation at the point.
  double PointValue() const;
};
} // namespace opensn
