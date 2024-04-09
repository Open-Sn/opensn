// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "framework/runtime.h"
#include "framework/physics/solver_base/solver.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/logging/log.h"
#include "field_functions_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(GetFieldFunctionHandleByName, fieldfunc, GetHandleByName);

int
GetFieldFunctionHandleByName(lua_State* L)
{
  const std::string fname = "fieldfunc.GetHandleByName";
  LuaCheckArgs<std::string>(L, fname);

  const auto ff_name = LuaArg<std::string>(L, 1);

  size_t ff_handle_counter = 0;
  std::vector<size_t> handles_that_matched;
  for (const auto& pff : opensn::field_function_stack)
  {
    if (pff->TextName() == ff_name)
      handles_that_matched.emplace_back(ff_handle_counter);
    ++ff_handle_counter;
  }

  size_t num_handles = handles_that_matched.size();

  if (num_handles == 0)
  {
    opensn::log.Log0Warning() << fname << ": No field-functions were found that "
                              << "matched the requested name:\"" << ff_name
                              << "\". A null handle will "
                              << "be returned." << std::endl;

    return LuaReturn(L);
  }

  if (num_handles > 1)
    opensn::log.Log0Warning() << fname << ": A total of " << num_handles
                              << " field-functions were found that matched the "
                              << " requested name. Only the first match will be "
                              << " returned.";

  auto handle = handles_that_matched.front();
  return LuaReturn(L, handle);
}

} // namespace opensnlua
