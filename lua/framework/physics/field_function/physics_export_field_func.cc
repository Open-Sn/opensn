#include "framework/lua.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "field_functions_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(ExportFieldFunctionToVTK, fieldfunc, ExportToVTK);
RegisterLuaFunctionNamespace(ExportMultiFieldFunctionToVTK, fieldfunc, ExportToVTKMulti);

int
ExportFieldFunctionToVTK(lua_State* L)
{
  const std::string fname = "fieldfunc.ExportToVTK";
  LuaCheckArgs<size_t, std::string>(L, fname);

  auto ff_handle = LuaArg<size_t>(L, 1);
  auto base_name = LuaArg<std::string>(L, 2);

  auto ff_base = opensn::GetStackItemPtr(opensn::field_function_stack, ff_handle, fname);
  auto ff = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff_base);

  OpenSnLogicalErrorIf(not ff, "Only grid-based field functions can be exported");

  //  ff->ExportToVTK(base_name);
  FieldFunctionGridBased::ExportMultipleToVTK(base_name, {ff});

  return 0;
}

int
ExportMultiFieldFunctionToVTK(lua_State* L)
{
  const std::string fname = "ExportMultiFieldFunctionToVTK";
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  auto base_name = LuaArg<std::string>(L, 2);

  LuaCheckTableValue(fname, L, 1);

  auto& ff_stack = opensn::field_function_stack;

  const size_t table_size = lua_rawlen(L, 1);
  std::vector<std::shared_ptr<const FieldFunctionGridBased>> ffs;
  ffs.reserve(table_size);
  for (int i = 0; i < table_size; ++i)
  {
    LuaPush(L, i + 1);
    lua_gettable(L, 1);

    std::shared_ptr<FieldFunction> ff_base = nullptr;
    if (lua_isinteger(L, -1))
    {
      int ff_handle = lua_tointeger(L, -1);
      lua_pop(L, 1);

      ff_base = opensn::GetStackItemPtr(ff_stack, ff_handle, fname);
    }
    else if (lua_isstring(L, -1))
    {
      const std::string ff_name = lua_tostring(L, -1);
      lua_pop(L, 1);

      for (auto& ff_ptr : ff_stack)
        if (ff_ptr->TextName() == ff_name)
        {
          ff_base = ff_ptr;
          break;
        }

      OpenSnInvalidArgumentIf(not ff_base,
                              "Field function with name \"" + ff_name + "\" could not be found.");
    }
    else
      OpenSnInvalidArgument("The field function specification can only be "
                            "string names or integer handles.");

    auto ff = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff_base);

    OpenSnLogicalErrorIf(not ff, "Only grid-based field functions can be exported");

    ffs.push_back(ff);
  } // for i

  FieldFunctionGridBased::ExportMultipleToVTK(base_name, ffs);

  return 0;
}

} // namespace opensnlua
