#ifndef CHI_FFINTER_LINE_H
#define CHI_FFINTER_LINE_H

#include "../chi_ffinterpolation.h"
#include "mesh/chi_mesh.h"

#include <petscksp.h>

namespace chi_mesh
{
  struct FieldFunctionContext
  {
    std::shared_ptr<chi_physics::FieldFunctionGridBased>    ref_ff;
    std::vector<double>            interpolation_points_values;
    std::vector<uint64_t>          interpolation_points_ass_cell;
    std::vector<bool>              interpolation_points_has_ass_cell;
  };
}

namespace chi_mesh
{
//###################################################################
/** A line based interpolation function.*/
class FieldFunctionInterpolationLine :
  public FieldFunctionInterpolation
{
protected:
  int                number_of_points_ = 2;
  chi_mesh::Vector3  pi_, pf_;

  std::vector<std::vector<double>>   custom_arrays_;
  std::vector<chi_mesh::Vector3>     interpolation_points_;
  std::vector<FieldFunctionContext>  ff_contexts_;

  double                         delta_d_ = 1.0;

public:
  FieldFunctionInterpolationLine() :
    FieldFunctionInterpolation(ff_interpolation::Type::LINE)
  {  }

  int& GetNumberOfPoints() {return number_of_points_;}
  chi_mesh::Vector3& GetInitialPoint() {return pi_;}
  chi_mesh::Vector3& GetFinalPoint() {return pf_;}
  std::vector<std::vector<double>>& GetCustomArrays() {return custom_arrays_;}
  std::vector<chi_mesh::Vector3>& GetInterpolationPoints()
  {return interpolation_points_;}
  std::vector<FieldFunctionContext>& GetFFContexts()
  {return ff_contexts_;}
  //01
  void Initialize() override;

  //02
  void Execute() override;

public:
  std::string GetDefaultFileBaseName() const override
  {return "ZLFFI";}
  void ExportPython(std::string base_name) override;
};
}//namespace chi_mesh



#endif