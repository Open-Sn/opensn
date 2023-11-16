#pragma once

#include "framework/mesh/field_function_interpolation/ffinterpolation.h"
#include "framework/mesh/logical_volume/logical_volume.h"

#include <petscksp.h>

namespace opensn
{

/**
 * Volume-wise field function interpolation.
 *
 * This interpolator allows the user to obtain quantities by logical
 * volume. If no logical volume is assigned to the method it will
 * default to operating over the entire volume.\n
 * \n
 * The method also supports a few primitive operations:
 *  - OP_VOLUME_AVG. Obtains the volume average of the field function
 *    of interest.
 *  - OP_VOLUME_SUM. Obtains the volume integral of the field function
 *    of interest.
 */
class FieldFunctionInterpolationVolume : public FieldFunctionInterpolation
{
protected:
  std::shared_ptr<LogicalVolume> logical_volume_ = nullptr;
  FieldFunctionInterpolationOperation op_type_ = FieldFunctionInterpolationOperation::OP_SUM;
  std::string op_lua_func_;
  double op_value_ = 0.0;

private:
  std::vector<uint64_t> cell_local_ids_inside_logvol_;

public:
  FieldFunctionInterpolationVolume()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::VOLUME)
  {
  }
  std::shared_ptr<LogicalVolume>& GetLogicalVolume() { return logical_volume_; }

  FieldFunctionInterpolationOperation& GetOperationType() { return op_type_; }

  std::string& GetOperationLuaFunction() { return op_lua_func_; }

  double& GetOpValue() { return op_value_; }

  void Initialize() override;
  void Execute() override;

#ifdef OPENSN_WITH_LUA
  /**
   * Calls the designated lua function
   */
  double CallLuaFunction(double ff_value, int mat_id) const;
#endif

  std::string GetDefaultFileBaseName() const override
  {
    return "ZVFFI";
  }
  void ExportPython(std::string base_name) override
  {
  }
};

} // namespace opensn
