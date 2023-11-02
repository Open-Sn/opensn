#ifndef CHI_FFINTER_SLICE_H
#define CHI_FFINTER_SLICE_H

#include "../chi_ffinterpolation.h"

#include "mesh/chi_mesh.h"

namespace chi_mesh
{
struct FFIFaceEdgeIntersection
{
  chi_mesh::Vector3 point;
  chi_mesh::Vector3 point2d;
  double point_value = 0.0;
};

struct FFICellIntersection
{
  uint64_t ref_cell_local_id = 0;
  std::vector<FFIFaceEdgeIntersection> intersections;
  chi_mesh::Vector3 intersection_centre;
  chi_mesh::Vector3 intersection_2d_centre;
  double cell_avg_value = 0.0;
};

}//namespace chi_mesh

//###################################################################
/** A slice based interpolation function.
 *
 * This functionality needs to cater for numerous spatial discretizations.
 * The most simple one is cell-averaged values and the more complicated ones
 * are PWLD and then CFEM.
 *
 * Cell average values requires computing the slice of the polyhedron and then
 * computing the centroid of that cut. This can be done cell by cell.*/
class chi_mesh::FieldFunctionInterpolationSlice : public chi_mesh::FieldFunctionInterpolation
{
protected:
  chi_mesh::Normal normal_ =chi_mesh::Normal(0.0, 0.0, 1.0);
  chi_mesh::Normal binorm_ =chi_mesh::Normal(0.0, 1.0, 0.0);
  chi_mesh::Normal tangent_=chi_mesh::Normal(1.0, 0.0, 0.0);
  chi_mesh::Vector3 plane_point_;

private:
  std::vector<FFICellIntersection>    cell_intersections_;
public:
  FieldFunctionInterpolationSlice() :
    FieldFunctionInterpolation(ff_interpolation::Type::SLICE)
  {  }

  chi_mesh::Normal& GetNormal() {return normal_;}
  chi_mesh::Normal& GetBiNorm() {return binorm_;}
  chi_mesh::Normal& GetTangent() {return tangent_;}
  chi_mesh::Vector3& GetPlanePoint() {return plane_point_;}
  //01
  void Initialize() override;

  //02
  void Execute() override;
public:
  //03
  std::string GetDefaultFileBaseName() const override
  {return "ZPFFI";}
  void ExportPython(std::string base_name) override;
};




#endif