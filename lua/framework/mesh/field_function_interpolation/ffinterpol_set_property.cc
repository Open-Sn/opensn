#include "framework/lua.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/physics/field_function/field_function_grid_based.h"

#include "framework/mesh/field_function_interpolation/ffinter_point.h"
#include "framework/mesh/field_function_interpolation/ffinter_slice.h"
#include "framework/mesh/field_function_interpolation/ffinter_line.h"
#include "framework/mesh/field_function_interpolation/ffinter_volume.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#define dcastPoint(x) dynamic_cast<FieldFunctionInterpolationPoint&>(x)
#define dcastLine(x) dynamic_cast<FieldFunctionInterpolationLine&>(x)
#define dcastSlice(x) dynamic_cast<FieldFunctionInterpolationSlice&>(x)
#define dcastVolume(x) dynamic_cast<FieldFunctionInterpolationVolume&>(x)

#include "ffinterpol_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiFFInterpolationSetProperty);
RegisterLuaConstantAsIs(PROBEPOINT, Varying(0));
RegisterLuaConstantAsIs(SLICE_POINT, Varying(1));
RegisterLuaConstantAsIs(SLICE_NORMAL, Varying(2));
RegisterLuaConstantAsIs(SLICE_TANGENT, Varying(3));
RegisterLuaConstantAsIs(SLICE_BINORM, Varying(4));
RegisterLuaConstantAsIs(OPERATION, Varying(5));
RegisterLuaConstantAsIs(OP_SUM, Varying(10));
RegisterLuaConstantAsIs(OP_AVG, Varying(11));
RegisterLuaConstantAsIs(OP_MAX, Varying(12));
RegisterLuaConstantAsIs(OP_SUM_LUA, Varying(13));
RegisterLuaConstantAsIs(OP_AVG_LUA, Varying(14));
RegisterLuaConstantAsIs(OP_MAX_LUA, Varying(15));
RegisterLuaConstantAsIs(LOGICAL_VOLUME, Varying(8));

RegisterLuaConstantAsIs(ADD_FIELDFUNCTION, Varying(9));
RegisterLuaConstantAsIs(SET_FIELDFUNCTIONS, Varying(10));

RegisterLuaConstantAsIs(LINE_FIRSTPOINT, Varying(11));
RegisterLuaConstantAsIs(LINE_SECONDPOINT, Varying(12));
RegisterLuaConstantAsIs(LINE_NUMBEROFPOINTS, Varying(13));
RegisterLuaConstantAsIs(LINE_CUSTOM_ARRAY, Varying(14));

int
chiFFInterpolationSetProperty(lua_State* L)
{
  const std::string fname = "chiFFInterpolationSetProperty";
  int numArgs = lua_gettop(L);

  // Get handle to field function
  const size_t ffihandle = lua_tonumber(L, 1);

  auto p_ffi =
    opensn::Chi::GetStackItemPtr(opensn::Chi::field_func_interpolation_stack, ffihandle, fname);

  // Process properties
  auto property = static_cast<FieldFunctionInterpolationProperty>(lua_tonumber(L, 2));
  // Check point properties
  if (property == FieldFunctionInterpolationProperty::PROBEPOINT)
    if (p_ffi->Type() != FieldFunctionInterpolationType::POINT)
      throw std::logic_error(
        "Point property" + std::to_string(static_cast<int>(property)) +
        " used in chiFFInterpolationSetProperty but FFI is not a point-probe.");

  // Check slice properties
  if ((property >= FieldFunctionInterpolationProperty::SLICEPOINT) &&
      (property <= FieldFunctionInterpolationProperty::SLICEBINORM))
    if (p_ffi->Type() != FieldFunctionInterpolationType::SLICE)
      throw std::logic_error("Slice property" + std::to_string(static_cast<int>(property)) +
                             " used in chiFFInterpolationSetProperty but FFI is not a slice.");

  // Check Line properties
  if ((property >= FieldFunctionInterpolationProperty::FIRSTPOINT) &&
      (property <= FieldFunctionInterpolationProperty::NUMBEROFPOINTS))
    if (p_ffi->Type() != FieldFunctionInterpolationType::LINE)
      throw std::logic_error("Line property " + std::to_string(static_cast<int>(property)) +
                             " used in chiFFInterpolationSetProperty but FFI is not a line.");

  // Generic
  if (property == FieldFunctionInterpolationProperty::ADD_FIELD_FUNCTION)
  {
    int ffhandle = lua_tonumber(L, 3);
    auto cur_ff_base =
      opensn::Chi::GetStackItemPtr(opensn::Chi::field_function_stack, ffhandle, fname);
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
      auto cur_ff_base =
        opensn::Chi::GetStackItemPtr(opensn::Chi::field_function_stack, ffhandle, fname);
      auto cur_ff = std::dynamic_pointer_cast<FieldFunctionGridBased>(cur_ff_base);

      p_ffi->GetFieldFunctions().push_back(cur_ff);
    } // for handle
  }
  else if (property == FieldFunctionInterpolationProperty::PROBEPOINT)
  {
    auto& cur_ffi = dcastPoint(*p_ffi);

    double x = lua_tonumber(L, 3);
    double y = lua_tonumber(L, 4);
    double z = lua_tonumber(L, 5);

    cur_ffi.GetPointOfInterest() = Vector3(x, y, z);
  }
  else if (property == FieldFunctionInterpolationProperty::SLICEPOINT)
  {
    auto& cur_ffi_slice = dcastSlice(*p_ffi);

    double x = lua_tonumber(L, 3);
    double y = lua_tonumber(L, 4);
    double z = lua_tonumber(L, 5);

    cur_ffi_slice.GetPlanePoint() = Vector3(x, y, z);
  }
  else if (property == FieldFunctionInterpolationProperty::SLICENORMAL)
  {
    auto& cur_ffi_slice = dcastSlice(*p_ffi);

    double x = lua_tonumber(L, 3);
    double y = lua_tonumber(L, 4);
    double z = lua_tonumber(L, 5);

    cur_ffi_slice.GetNormal() = Vector3(x, y, z).Normalized();
  }
  else if (property == FieldFunctionInterpolationProperty::SLICETANGENT)
  {
    auto& cur_ffi_slice = dcastSlice(*p_ffi);

    double x = lua_tonumber(L, 3);
    double y = lua_tonumber(L, 4);
    double z = lua_tonumber(L, 5);

    cur_ffi_slice.GetTangent() = Vector3(x, y, z).Normalized();
  }
  else if (property == FieldFunctionInterpolationProperty::SLICEBINORM)
  {
    auto& cur_ffi_slice = dcastSlice(*p_ffi);

    double x = lua_tonumber(L, 3);
    double y = lua_tonumber(L, 4);
    double z = lua_tonumber(L, 5);

    cur_ffi_slice.GetBiNorm() = Vector3(x, y, z).Normalized();
  }
  else if (property == FieldFunctionInterpolationProperty::FIRSTPOINT)
  {
    if (numArgs != 5) LuaPostArgAmountError("chiFFInterpolationSetProperty", 5, numArgs);

    auto& cur_ffi_line = dcastLine(*p_ffi);

    Vector3 point(lua_tonumber(L, 3), lua_tonumber(L, 4), lua_tonumber(L, 5));
    cur_ffi_line.GetInitialPoint() = point;
  }
  else if (property == FieldFunctionInterpolationProperty::SECONDPOINT)
  {
    if (numArgs != 5) LuaPostArgAmountError("chiFFInterpolationSetProperty", 5, numArgs);

    auto& cur_ffi_line = dcastLine(*p_ffi);

    Vector3 point(lua_tonumber(L, 3), lua_tonumber(L, 4), lua_tonumber(L, 5));
    cur_ffi_line.GetFinalPoint() = point;
  }
  else if (property == FieldFunctionInterpolationProperty::NUMBEROFPOINTS)
  {
    if (numArgs != 3) LuaPostArgAmountError("chiFFInterpolationSetProperty", 3, numArgs);

    auto& cur_ffi_line = dcastLine(*p_ffi);

    int num_points = lua_tonumber(L, 3);

    if (num_points < 2)
    {
      opensn::log.LogAllError() << "Line property FFI_LINE_NUMBEROFPOINTS"
                                << " used in chiFFInterpolationSetProperty. Number of points must"
                                << " be greater than or equal to 2.";
      opensn::Chi::Exit(EXIT_FAILURE);
    }
    cur_ffi_line.GetNumberOfPoints() = num_points;
  }
  else if (property == FieldFunctionInterpolationProperty::CUSTOM_ARRAY)
  {
    if (numArgs != 3) LuaPostArgAmountError("chiFFInterpolationSetProperty", 3, numArgs);

    auto& cur_ffi_line = dcastLine(*p_ffi);

    if (not lua_istable(L, 3))
    {
      opensn::log.LogAllError() << "Line property FFI_LINE_CUSTOM_ARRAY"
                                << " used in chiFFInterpolationSetProperty. Argument 3 is expected "
                                   "to be an array.";
      opensn::Chi::Exit(EXIT_FAILURE);
    }

    const size_t table_len = lua_rawlen(L, 3);

    std::vector<double> new_array(table_len, 0.0);
    for (int k = 0; k < table_len; ++k)
    {
      lua_pushnumber(L, k + 1);
      lua_gettable(L, 3);
      new_array[k] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }

    cur_ffi_line.GetCustomArrays().push_back(new_array);
  }
  else if (property == FieldFunctionInterpolationProperty::OPERATION)
  {
    if (numArgs != 3 and numArgs != 4)
      LuaPostArgAmountError("chiFFInterpolationSetProperty", 3, numArgs);

    if (p_ffi->Type() != FieldFunctionInterpolationType::VOLUME)
      throw std::logic_error("Volume property FFI_PROP_OPERATION"
                             " used in chiFFInterpolationSetProperty can only be used with "
                             "Volume type interpolations.");

    auto& cur_ffi_volume = dcastVolume(*p_ffi);

    int op_type = lua_tonumber(L, 3);

    int OP_SUM = static_cast<int>(FieldFunctionInterpolationOperation::OP_SUM);
    int OP_MAX_LUA = static_cast<int>(FieldFunctionInterpolationOperation::OP_MAX_LUA);
    int OP_SUM_LUA = static_cast<int>(FieldFunctionInterpolationOperation::OP_SUM_LUA);

    if (!((op_type >= OP_SUM) && (op_type <= OP_MAX_LUA)))
    {
      opensn::log.LogAllError() << "Volume property FFI_PROP_OPERATION"
                                << " used in chiFFInterpolationSetProperty. Unsupported OPERATON."
                                << " Supported types are OP_AVG and OP_SUM. " << op_type;
      opensn::Chi::Exit(EXIT_FAILURE);
    }

    if ((op_type >= OP_SUM_LUA) and (op_type <= OP_MAX_LUA))
    {
      if (numArgs != 4) LuaPostArgAmountError("chiFFInterpolationSetProperty", 4, numArgs);

      const char* func_name = lua_tostring(L, 4);
      cur_ffi_volume.GetOperationLuaFunction() = std::string(func_name);
    }

    cur_ffi_volume.GetOperationType() = static_cast<FieldFunctionInterpolationOperation>(op_type);
  }
  else if (property == FieldFunctionInterpolationProperty::LOGICAL_VOLUME)
  {
    if (numArgs != 3) LuaPostArgAmountError("chiFFInterpolationSetProperty", 3, numArgs);

    int logvol_hndle = lua_tonumber(L, 3);

    auto p_logical_volume = std::dynamic_pointer_cast<LogicalVolume>(
      opensn::Chi::GetStackItemPtr(opensn::Chi::object_stack, logvol_hndle, fname));

    if (p_ffi->Type() != FieldFunctionInterpolationType::VOLUME)
      throw std::logic_error("Volume property FFI_PROP_LOGICAL_VOLUME"
                             " used in chiFFInterpolationSetProperty can only be used with "
                             "Volume type interpolations.");

    auto& cur_ffi_volume = dcastVolume(*p_ffi);

    cur_ffi_volume.GetLogicalVolume() = p_logical_volume;
  }
  else // Fall back
  {
    opensn::log.LogAllError() << "Invalid PropertyIndex used in chiFFInterpolationSetProperty.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }

  return 0;
}
