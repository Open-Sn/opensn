#pragma once

#include "framework/mesh/FieldFunctionInterpolation/chi_ffinterpolation.h"
#include "framework/mesh/chi_mesh.h"

namespace chi_mesh
{

/** A line based interpolation function.*/
class FieldFunctionInterpolationPoint : public FieldFunctionInterpolation
{
protected:
  chi_mesh::Vector3 point_of_interest_;

  bool locally_owned_ = false;
  uint64_t owning_cell_gid_ = 0;
  double point_value_ = 0.0;

public:
  FieldFunctionInterpolationPoint() : FieldFunctionInterpolation(ff_interpolation::Type::POINT) {}

  chi_mesh::Vector3& GetPointOfInterest() { return point_of_interest_; }

  void Initialize() override;
  void Execute() override;
  /**Gets the value of the field function evaluation at the point.*/
  double GetPointValue() const;

public:
  std::string GetDefaultFileBaseName() const override { return ""; }
  void ExportPython(std::string base_name) override{};
};
} // namespace chi_mesh
