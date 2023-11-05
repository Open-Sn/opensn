#include "chi_lua.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "fieldfunctions_lua.h"
#include "console/chi_console.h"

RegisterLuaFunctionAsIs(chiExportFieldFunctionToVTK);
RegisterLuaFunctionAsIs(chiExportMultiFieldFunctionToVTK);

// #############################################################################
/** Exports a field function to VTK format.
 *
\param FFHandle int Global handle to the field function.
\param BaseName char Base name for the exported file.

\ingroup LuaFieldFunc
\author Jan*/
int
chiExportFieldFunctionToVTK(lua_State* L)
{
  const std::string fname = "chiExportFieldFunctionToVTK";
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  int ff_handle = lua_tonumber(L, 1);
  const char* base_name = lua_tostring(L, 2);

  typedef chi_physics::FieldFunctionGridBased FFGridBased;

  auto ff_base = Chi::GetStackItemPtr(Chi::field_function_stack, ff_handle, fname);
  auto ff = std::dynamic_pointer_cast<FFGridBased>(ff_base);

  ChiLogicalErrorIf(not ff, "Only grid-based field functions can be exported");

  //  ff->ExportToVTK(base_name);
  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK(base_name, {ff});

  return 0;
}

// #############################################################################
/** Exports all the field functions in a list to VTK format.
 *  *
\param listFFHandles table Global handles or names to the field functions
\param BaseName char Base name for the exported file.

\ingroup LuaFieldFunc
\author Jan*/
int
chiExportMultiFieldFunctionToVTK(lua_State* L)
{
  const std::string fname = "chiExportMultiFieldFunctionToVTK";
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  const char* base_name = lua_tostring(L, 2);

  LuaCheckTableValue(fname, L, 1);

  auto& ff_stack = Chi::field_function_stack;

  const size_t table_size = lua_rawlen(L, 1);
  std::vector<std::shared_ptr<const chi_physics::FieldFunctionGridBased>> ffs;
  ffs.reserve(table_size);
  for (int i = 0; i < table_size; ++i)
  {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, 1);

    std::shared_ptr<chi_physics::FieldFunction> ff_base = nullptr;
    if (lua_isinteger(L, -1))
    {
      int ff_handle = lua_tonumber(L, -1);
      lua_pop(L, 1);

      ff_base = Chi::GetStackItemPtr(ff_stack, ff_handle, fname);
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

      ChiInvalidArgumentIf(not ff_base,
                           "Field function with name \"" + ff_name + "\" could not be found.");
    }
    else
      ChiInvalidArgument("The field function specification can only be "
                         "string names or integer handles.");

    typedef chi_physics::FieldFunctionGridBased FFGridBased;
    auto ff = std::dynamic_pointer_cast<FFGridBased>(ff_base);

    ChiLogicalErrorIf(not ff, "Only grid-based field functions can be exported");

    ffs.push_back(ff);
  } // for i

  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK(base_name, ffs);

  return 0;
}
