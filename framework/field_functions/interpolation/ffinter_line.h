#pragma once

#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/mesh/mesh.h"

#include <petscksp.h>

namespace opensn
{

struct FieldFunctionContext
{
  std::shared_ptr<FieldFunctionGridBased> ref_ff;
  std::vector<double> interpolation_points_values;
  std::vector<uint64_t> interpolation_points_ass_cell;
  std::vector<bool> interpolation_points_has_ass_cell;
};

/** A line based interpolation function.*/
class FieldFunctionInterpolationLine : public FieldFunctionInterpolation
{
protected:
  int number_of_points_ = 2;
  Vector3 pi_, pf_;

  std::vector<std::vector<double>> custom_arrays_;
  std::vector<Vector3> interpolation_points_;
  std::vector<FieldFunctionContext> ff_contexts_;

  double delta_d_ = 1.0;

public:
  FieldFunctionInterpolationLine()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::LINE)
  {
  }

  virtual ~FieldFunctionInterpolationLine() {}

  int& GetNumberOfPoints() { return number_of_points_; }
  Vector3& GetInitialPoint() { return pi_; }
  Vector3& GetFinalPoint() { return pf_; }
  std::vector<std::vector<double>>& GetCustomArrays() { return custom_arrays_; }
  std::vector<Vector3>& GetInterpolationPoints() { return interpolation_points_; }
  std::vector<FieldFunctionContext>& GetFFContexts() { return ff_contexts_; }
  void Initialize() override;
  void Execute() override;

public:
  std::string GetDefaultFileBaseName() const override { return "ZLFFI"; }
  void ExportPython(std::string base_name) override;
};

} // namespace opensn
