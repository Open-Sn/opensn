#include "framework/lua.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinter_point.h"
#include "framework/field_functions/interpolation/ffinter_slice.h"
#include "framework/field_functions/interpolation/ffinter_line.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include "lua/framework/math/functions/lua_scalar_material_function.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "ffinterpol_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(FFInterpolationSetProperty, fieldfunc, SetProperty);
RegisterLuaConstantAsIs(PROBEPOINT, Varying(0));
RegisterLuaConstantAsIs(SLICE_POINT, Varying(1));
RegisterLuaConstantAsIs(SLICE_NORMAL, Varying(2));
RegisterLuaConstantAsIs(SLICE_TANGENT, Varying(3));
RegisterLuaConstantAsIs(SLICE_BINORM, Varying(4));
RegisterLuaConstantAsIs(OPERATION, Varying(5));
RegisterLuaConstantAsIs(OP_SUM, Varying(10));
RegisterLuaConstantAsIs(OP_AVG, Varying(11));
RegisterLuaConstantAsIs(OP_MAX, Varying(12));
RegisterLuaConstantAsIs(OP_SUM_FUNC, Varying(13));
RegisterLuaConstantAsIs(OP_AVG_FUNC, Varying(14));
RegisterLuaConstantAsIs(OP_MAX_FUNC, Varying(15));
RegisterLuaConstantAsIs(LOGICAL_VOLUME, Varying(8));

RegisterLuaConstantAsIs(ADD_FIELDFUNCTION, Varying(9));
RegisterLuaConstantAsIs(SET_FIELDFUNCTIONS, Varying(10));

RegisterLuaConstantAsIs(LINE_FIRSTPOINT, Varying(11));
RegisterLuaConstantAsIs(LINE_SECONDPOINT, Varying(12));
RegisterLuaConstantAsIs(LINE_NUMBEROFPOINTS, Varying(13));
RegisterLuaConstantAsIs(LINE_CUSTOM_ARRAY, Varying(14));

namespace
{

std::shared_ptr<LuaScalarMaterialFunction>
CreateFunction(const std::string& function_name)
{
  ParameterBlock blk;
  blk.AddParameter("lua_function_name", function_name);
  InputParameters params = LuaScalarMaterialFunction::GetInputParameters();
  params.AssignParameters(blk);
  return std::make_shared<LuaScalarMaterialFunction>(params);
}

} // namespace

int
FFInterpolationSetProperty(lua_State* L)
{
  const std::string fname = "FFInterpolationSetProperty";
  int numArgs = lua_gettop(L);

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);

  auto p_ffi = opensn::GetStackItemPtr(opensn::field_func_interpolation_stack, ffihandle, fname);

  // Process properties
  auto property = static_cast<FieldFunctionInterpolationProperty>(LuaArg<int>(L, 2));
  // Check point properties
  if (property == FieldFunctionInterpolationProperty::PROBEPOINT)
    if (p_ffi->Type() != FieldFunctionInterpolationType::POINT)
      throw std::logic_error("Point property" + std::to_string(static_cast<int>(property)) +
                             " used in FFInterpolationSetProperty but FFI is not a point-probe.");

  // Check slice properties
  if ((property >= FieldFunctionInterpolationProperty::SLICEPOINT) and
      (property <= FieldFunctionInterpolationProperty::SLICEBINORM))
    if (p_ffi->Type() != FieldFunctionInterpolationType::SLICE)
      throw std::logic_error("Slice property" + std::to_string(static_cast<int>(property)) +
                             " used in FFInterpolationSetProperty but FFI is not a slice.");

  // Check Line properties
  if ((property >= FieldFunctionInterpolationProperty::FIRSTPOINT) and
      (property <= FieldFunctionInterpolationProperty::NUMBEROFPOINTS))
    if (p_ffi->Type() != FieldFunctionInterpolationType::LINE)
      throw std::logic_error("Line property " + std::to_string(static_cast<int>(property)) +
                             " used in FFInterpolationSetProperty but FFI is not a line.");

  // Generic
  if (property == FieldFunctionInterpolationProperty::ADD_FIELD_FUNCTION)
  {
    auto ffhandle = LuaArg<size_t>(L, 3);
    auto cur_ff_base = opensn::GetStackItemPtr(opensn::field_function_stack, ffhandle, fname);
    auto cur_ff = std::dynamic_pointer_cast<FieldFunctionGridBased>(cur_ff_base);

    p_ffi->GetFieldFunctions().push_back(cur_ff);
  }
  else if (property == FieldFunctionInterpolationProperty::SET_FIELD_FUNCTIONS)
  {
    LuaCheckTableValue(fname, L, 3);
    std::vector<double> handle_array;
    LuaPopulateVectorFrom1DArray(fname, L, 3, handle_array);

    for (double handle_d : handle_array)
    {
      const auto ffhandle = static_cast<int>(handle_d);
      auto cur_ff_base = opensn::GetStackItemPtr(opensn::field_function_stack, ffhandle, fname);
      auto cur_ff = std::dynamic_pointer_cast<FieldFunctionGridBased>(cur_ff_base);

      p_ffi->GetFieldFunctions().push_back(cur_ff);
    } // for handle
  }
  else if (property == FieldFunctionInterpolationProperty::PROBEPOINT)
  {
    auto& cur_ffi = dynamic_cast<FieldFunctionInterpolationPoint&>(*p_ffi);
    cur_ffi.GetPointOfInterest() = LuaArg<Vector3>(L, 3);
  }
  else if (property == FieldFunctionInterpolationProperty::SLICEPOINT)
  {
    auto& cur_ffi_slice = dynamic_cast<FieldFunctionInterpolationSlice&>(*p_ffi);
    cur_ffi_slice.GetPlanePoint() = LuaArg<Vector3>(L, 3);
  }
  else if (property == FieldFunctionInterpolationProperty::SLICENORMAL)
  {
    auto& cur_ffi_slice = dynamic_cast<FieldFunctionInterpolationSlice&>(*p_ffi);
    cur_ffi_slice.GetNormal() = LuaArg<Vector3>(L, 3).Normalized();
  }
  else if (property == FieldFunctionInterpolationProperty::SLICETANGENT)
  {
    auto& cur_ffi_slice = dynamic_cast<FieldFunctionInterpolationSlice&>(*p_ffi);
    cur_ffi_slice.GetTangent() = LuaArg<Vector3>(L, 3).Normalized();
  }
  else if (property == FieldFunctionInterpolationProperty::SLICEBINORM)
  {
    auto& cur_ffi_slice = dynamic_cast<FieldFunctionInterpolationSlice&>(*p_ffi);
    cur_ffi_slice.GetBiNorm() = LuaArg<Vector3>(L, 3).Normalized();
  }
  else if (property == FieldFunctionInterpolationProperty::FIRSTPOINT)
  {
    if (numArgs != 3)
      LuaPostArgAmountError("FFInterpolationSetProperty", 3, numArgs);

    auto& cur_ffi_line = dynamic_cast<FieldFunctionInterpolationLine&>(*p_ffi);
    cur_ffi_line.GetInitialPoint() = LuaArg<Vector3>(L, 3);
  }
  else if (property == FieldFunctionInterpolationProperty::SECONDPOINT)
  {
    if (numArgs != 3)
      LuaPostArgAmountError("FFInterpolationSetProperty", 3, numArgs);

    auto& cur_ffi_line = dynamic_cast<FieldFunctionInterpolationLine&>(*p_ffi);
    cur_ffi_line.GetFinalPoint() = LuaArg<Vector3>(L, 3);
  }
  else if (property == FieldFunctionInterpolationProperty::NUMBEROFPOINTS)
  {
    if (numArgs != 3)
      LuaPostArgAmountError("FFInterpolationSetProperty", 3, numArgs);

    auto& cur_ffi_line = dynamic_cast<FieldFunctionInterpolationLine&>(*p_ffi);

    auto num_points = LuaArg<int>(L, 3);

    if (num_points < 2)
    {
      opensn::log.LogAllError() << "Line property FFI_LINE_NUMBEROFPOINTS"
                                << " used in FFInterpolationSetProperty. Number of points must"
                                << " be greater than or equal to 2.";
      opensn::Exit(EXIT_FAILURE);
    }
    cur_ffi_line.GetNumberOfPoints() = num_points;
  }
  else if (property == FieldFunctionInterpolationProperty::CUSTOM_ARRAY)
  {
    if (numArgs != 3)
      LuaPostArgAmountError("FFInterpolationSetProperty", 3, numArgs);

    auto& cur_ffi_line = dynamic_cast<FieldFunctionInterpolationLine&>(*p_ffi);
    auto new_array = LuaArg<std::vector<double>>(L, 3);
    cur_ffi_line.GetCustomArrays().push_back(new_array);
  }
  else if (property == FieldFunctionInterpolationProperty::OPERATION)
  {
    if (numArgs != 3 and numArgs != 4)
      LuaPostArgAmountError("FFInterpolationSetProperty", 3, numArgs);

    if (p_ffi->Type() != FieldFunctionInterpolationType::VOLUME)
      throw std::logic_error("Volume property FFI_PROP_OPERATION"
                             " used in FFInterpolationSetProperty can only be used with "
                             "Volume type interpolations.");

    auto& cur_ffi_volume = dynamic_cast<FieldFunctionInterpolationVolume&>(*p_ffi);

    auto op_type = LuaArg<int>(L, 3);

    int OP_SUM = static_cast<int>(FieldFunctionInterpolationOperation::OP_SUM);
    int OP_MAX_FUNC = static_cast<int>(FieldFunctionInterpolationOperation::OP_MAX_FUNC);
    int OP_SUM_FUNC = static_cast<int>(FieldFunctionInterpolationOperation::OP_SUM_FUNC);

    if (not((op_type >= OP_SUM) and (op_type <= OP_MAX_FUNC)))
    {
      opensn::log.LogAllError() << "Volume property FFI_PROP_OPERATION"
                                << " used in FFInterpolationSetProperty. Unsupported OPERATON."
                                << " Supported types are OP_AVG and OP_SUM. " << op_type;
      opensn::Exit(EXIT_FAILURE);
    }

    if ((op_type >= OP_SUM_FUNC) and (op_type <= OP_MAX_FUNC))
    {
      if (numArgs != 4)
        LuaPostArgAmountError("FFInterpolationSetProperty", 4, numArgs);
      const auto func_name = LuaArg<std::string>(L, 4);
      auto operation_function = CreateFunction(func_name);
      opensn::function_stack.push_back(operation_function);
      cur_ffi_volume.SetOperationFunction(operation_function);
    }

    cur_ffi_volume.GetOperationType() = static_cast<FieldFunctionInterpolationOperation>(op_type);
  }
  else if (property == FieldFunctionInterpolationProperty::LOGICAL_VOLUME)
  {
    if (numArgs != 3)
      LuaPostArgAmountError("FFInterpolationSetProperty", 3, numArgs);

    auto logvol_hndle = LuaArg<int>(L, 3);

    auto p_logical_volume = std::dynamic_pointer_cast<LogicalVolume>(
      opensn::GetStackItemPtr(opensn::object_stack, logvol_hndle, fname));

    if (p_ffi->Type() != FieldFunctionInterpolationType::VOLUME)
      throw std::logic_error("Volume property FFI_PROP_LOGICAL_VOLUME"
                             " used in FFInterpolationSetProperty can only be used with "
                             "Volume type interpolations.");

    auto& cur_ffi_volume = dynamic_cast<FieldFunctionInterpolationVolume&>(*p_ffi);

    cur_ffi_volume.GetLogicalVolume() = p_logical_volume;
  }
  else // Fall back
  {
    opensn::log.LogAllError() << "Invalid PropertyIndex used in FFInterpolationSetProperty.";
    opensn::Exit(EXIT_FAILURE);
  }

  return 0;
}

} // namespace opensnlua
