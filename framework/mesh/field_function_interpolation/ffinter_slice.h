#pragma once

#include "framework/mesh/field_function_interpolation/ffinterpolation.h"
#include "framework/mesh/mesh.h"

namespace opensn
{

struct FFIFaceEdgeIntersection
{
  Vector3 point;
  Vector3 point2d;
  double point_value = 0.0;
};

struct FFICellIntersection
{
  uint64_t ref_cell_local_id = 0;
  std::vector<FFIFaceEdgeIntersection> intersections;
  Vector3 intersection_centre;
  Vector3 intersection_2d_centre;
  double cell_avg_value = 0.0;
};

/**
 * A slice based interpolation function.
 *
 * This functionality needs to cater for numerous spatial discretizations.
 * The most simple one is cell-averaged values and the more complicated ones
 * are PWLD and then CFEM.
 *
 * Cell average values requires computing the slice of the polyhedron and then
 * computing the centroid of that cut. This can be done cell by cell.
 */
class FieldFunctionInterpolationSlice : public FieldFunctionInterpolation
{
protected:
  Normal normal_ = Normal(0.0, 0.0, 1.0);
  Normal binorm_ = Normal(0.0, 1.0, 0.0);
  Normal tangent_ = Normal(1.0, 0.0, 0.0);
  Vector3 plane_point_;

private:
  std::vector<FFICellIntersection> cell_intersections_;

public:
  FieldFunctionInterpolationSlice()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::SLICE)
  {
  }

  Normal& GetNormal() { return normal_; }
  Normal& GetBiNorm() { return binorm_; }
  Normal& GetTangent() { return tangent_; }
  Vector3& GetPlanePoint() { return plane_point_; }

  /**
   * Initializes the data structures necessary for interpolation. This is independent of the physics
   * and hence is a routine on its own.
   *
   * The first step of this initialization is to determine which cells are intersected by this
   * plane. For polyhedrons this is evaluated tet-by-tet.
   *
   * The second step is find where face-edges are intersected. This will effectively create
   * intersection polygons.
   */
  void Initialize() override;
  void Execute() override;

  std::string GetDefaultFileBaseName() const override { return "ZPFFI"; }
  void ExportPython(std::string base_name) override;
};

} // namespace opensn
