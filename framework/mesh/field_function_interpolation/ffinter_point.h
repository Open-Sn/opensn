#pragma once

#include "framework/mesh/field_function_interpolation/ffinterpolation.h"
#include "framework/mesh/mesh.h"

namespace opensn
{

/** A line based interpolation function.*/
class FieldFunctionInterpolationPoint : public FieldFunctionInterpolation
{
protected:
  Vector3 point_of_interest_;

  bool locally_owned_ = false;
  uint64_t owning_cell_gid_ = 0;
  double point_value_ = 0.0;

public:
  FieldFunctionInterpolationPoint()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::POINT)
  {
  }

  Vector3& GetPointOfInterest() { return point_of_interest_; }

  void Initialize() override;
  void Execute() override;
  /**Gets the value of the field function evaluation at the point.*/
  double GetPointValue() const;

public:
  std::string GetDefaultFileBaseName() const override { return ""; }
  void ExportPython(std::string base_name) override{};
};
} // namespace opensn
