// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinter_point.h"
#include "framework/field_functions/interpolation/ffinter_slice.h"
#include "framework/field_functions/interpolation/ffinter_line.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include "lua/framework/math/functions/lua_scalar_material_function.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/framework/mesh/field_function_interpolation/ffinterpol.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(FFInterpolationSetProperty, fieldfunc, SetProperty);
RegisterLuaConstant(PROBEPOINT, Varying(0));
RegisterLuaConstant(SLICE_POINT, Varying(1));
RegisterLuaConstant(SLICE_NORMAL, Varying(2));
RegisterLuaConstant(SLICE_TANGENT, Varying(3));
RegisterLuaConstant(SLICE_BINORM, Varying(4));
RegisterLuaConstant(OPERATION, Varying(5));
RegisterLuaConstant(OP_SUM, Varying(10));
RegisterLuaConstant(OP_AVG, Varying(11));
RegisterLuaConstant(OP_MAX, Varying(12));
RegisterLuaConstant(OP_SUM_FUNC, Varying(13));
RegisterLuaConstant(OP_AVG_FUNC, Varying(14));
RegisterLuaConstant(OP_MAX_FUNC, Varying(15));
RegisterLuaConstant(LOGICAL_VOLUME, Varying(8));
RegisterLuaConstant(ADD_FIELDFUNCTION, Varying(9));
RegisterLuaConstant(SET_FIELDFUNCTIONS, Varying(10));
RegisterLuaConstant(LINE_FIRSTPOINT, Varying(11));
RegisterLuaConstant(LINE_SECONDPOINT, Varying(12));
RegisterLuaConstant(LINE_NUMBEROFPOINTS, Varying(13));
RegisterLuaConstant(LINE_CUSTOM_ARRAY, Varying(14));

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
  const std::string fname = "fieldfunc.SetProperty";
  LuaCheckArgs<size_t, int>(L, fname);

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
    auto handle_array = LuaArg<std::vector<double>>(L, 3);
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
    LuaCheckArgs<size_t, int, Vector3>(L, fname);

    auto& cur_ffi_line = dynamic_cast<FieldFunctionInterpolationLine&>(*p_ffi);
    cur_ffi_line.SetInitialPoint(LuaArg<Vector3>(L, 3));
  }
  else if (property == FieldFunctionInterpolationProperty::SECONDPOINT)
  {
    LuaCheckArgs<size_t, int, Vector3>(L, fname);

    auto& cur_ffi_line = dynamic_cast<FieldFunctionInterpolationLine&>(*p_ffi);
    cur_ffi_line.SetFinalPoint(LuaArg<Vector3>(L, 3));
  }
  else if (property == FieldFunctionInterpolationProperty::NUMBEROFPOINTS)
  {
    LuaCheckArgs<size_t, int, int>(L, fname);

    auto& cur_ffi_line = dynamic_cast<FieldFunctionInterpolationLine&>(*p_ffi);
    auto num_points = LuaArg<int>(L, 3);
    if (num_points < 2)
    {
      throw std::logic_error("Line property FFI_LINE_NUMBEROFPOINTS used in "
                             "FFInterpolationSetProperty. Number of points must be greater than or "
                             "equal to 2.");
    }
    cur_ffi_line.SetNumberOfPoints(num_points);
  }
  else if (property == FieldFunctionInterpolationProperty::OPERATION)
  {
    LuaCheckArgs<size_t, int, int>(L, fname);
    auto op_type = LuaArg<int>(L, 3);

    int OP_SUM = static_cast<int>(FieldFunctionInterpolationOperation::OP_SUM);
    int OP_AVG = static_cast<int>(FieldFunctionInterpolationOperation::OP_AVG);
    int OP_MAX = static_cast<int>(FieldFunctionInterpolationOperation::OP_MAX);
    int OP_SUM_FUNC = static_cast<int>(FieldFunctionInterpolationOperation::OP_SUM_FUNC);
    int OP_AVG_FUNC = static_cast<int>(FieldFunctionInterpolationOperation::OP_AVG_FUNC);
    int OP_MAX_FUNC = static_cast<int>(FieldFunctionInterpolationOperation::OP_MAX_FUNC);

    if (p_ffi->Type() == FieldFunctionInterpolationType::VOLUME)
    {
      if (op_type < OP_SUM or op_type > OP_MAX_FUNC)
      {
        throw std::logic_error("FFI_PROP_OPERATION used in FFInterpolationSetProperty. Unsupported "
                               "operation type. Supported types are OP_SUM, OP_AVG, OP_MAX, "
                               "OP_SUM_FUNC, OP_AVG_FUNC and OP_MAX_FUNC.");
      }

      if (op_type >= OP_SUM_FUNC)
      {
        LuaCheckArgs<size_t, int, int, std::string>(L, fname);
        const auto func_name = LuaArg<std::string>(L, 4);
        auto operation_function = CreateFunction(func_name);
        opensn::function_stack.push_back(operation_function);
        auto& cur_ffi_volume = dynamic_cast<FieldFunctionInterpolationVolume&>(*p_ffi);
        cur_ffi_volume.SetOperationFunction(operation_function);
      }

      auto& cur_ffi_volume = dynamic_cast<FieldFunctionInterpolationVolume&>(*p_ffi);
      cur_ffi_volume.GetOperationType() = static_cast<FieldFunctionInterpolationOperation>(op_type);
    }
    else if (p_ffi->Type() == FieldFunctionInterpolationType::LINE)
    {
      if (op_type < OP_SUM && op_type > OP_MAX)
      {
        throw std::logic_error("FFI_PROP_OPERATION used in FFInterpolationSetProperty. Unsupported "
                               "operation type. Supported types are OP_SUM, OP_AVG, or OP_MAX.");
      }

      auto& cur_ffi_line = dynamic_cast<FieldFunctionInterpolationLine&>(*p_ffi);
      cur_ffi_line.SetOperationType(static_cast<FieldFunctionInterpolationOperation>(op_type));
    }
    else
    {
      throw std::logic_error("FFI_PROP_OPERATION used in FFInterpolationSetProperty can only be "
                             "used with volume or line interpolations.");
    }
  }
  else if (property == FieldFunctionInterpolationProperty::LOGICAL_VOLUME)
  {
    LuaCheckArgs<size_t, int, int>(L, fname);

    auto logvol_hndle = LuaArg<int>(L, 3);
    auto p_logical_volume = std::dynamic_pointer_cast<LogicalVolume>(
      opensn::GetStackItemPtr(opensn::object_stack, logvol_hndle, fname));

    if (p_ffi->Type() != FieldFunctionInterpolationType::VOLUME)
    {
      throw std::logic_error("Volume property FFI_PROP_LOGICAL_VOLUME use in "
                             "FFInterpolationSetProperty can only be used with volume type "
                             "interpolations.");
    }

    auto& cur_ffi_volume = dynamic_cast<FieldFunctionInterpolationVolume&>(*p_ffi);
    cur_ffi_volume.GetLogicalVolume() = p_logical_volume;
  }
  else
    throw std::logic_error("Invalid PropertyIndex used in FFInterpolationSetProperty.");

  return LuaReturn(L);
}

} // namespace opensnlua
