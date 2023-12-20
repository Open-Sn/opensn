#include "physics_solver_lua_utils.h"

#include "framework/physics/solver_base/solver.h"
#include "framework/physics/field_function/field_function_grid_based.h"
#include "framework/physics/physics_event_publisher.h"

#include "framework/object_factory.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{
RegisterLuaFunctionAsIs(chiSolverCreate);

RegisterLuaFunctionAsIs(chiSolverInitialize);
RegisterLuaFunctionAsIs(chiSolverExecute);
RegisterLuaFunctionAsIs(chiSolverStep);
RegisterLuaFunctionAsIs(chiSolverAdvance);
RegisterLuaFunctionAsIs(chiSolverSetBasicOption);
RegisterLuaFunctionAsIs(chiSolverGetName);
RegisterLuaFunctionAsIs(chiSolverGetFieldFunctionList);
RegisterLuaFunctionAsIs(chiSolverGetInfo);
RegisterLuaFunctionAsIs(chiSolverSetProperties);

int
chiSolverCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckTableValue(fname, L, 1);

  const auto params = TableParserAsParameterBlock::ParseTable(L, 1);

  const auto& object_maker = ChiObjectFactory::GetInstance();
  const size_t handle = object_maker.MakeRegisteredObject(params);

  lua_pushinteger(L, static_cast<lua_Integer>(handle));
  return 1;
}

int
chiSolverInitialize(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  PhysicsEventPublisher::GetInstance().SolverInitialize(solver);

  return 0;
}

int
chiSolverExecute(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  PhysicsEventPublisher::GetInstance().SolverExecute(solver);

  return 0;
}

int
chiSolverStep(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  PhysicsEventPublisher::GetInstance().SolverStep(solver);

  return 0;
}

int
chiSolverAdvance(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  PhysicsEventPublisher::GetInstance().SolverAdvance(solver);

  return 0;
}

int
chiSolverSetBasicOption(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const int solver_handle = lua_tointeger(L, 1);
  const std::string option_name = lua_tostring(L, 2);

  auto& solver = opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  try
  {
    auto& option = solver.GetBasicOptions()[option_name];

    switch (option.Type())
    {
      case VaryingDataType::VOID:
      case VaryingDataType::ARBITRARY_BYTES:
        throw std::logic_error("Solver:" + solver.TextName() + " option:" + option_name +
                               " is of invalid type."
                               " This indicates an implementation problem.");
      case VaryingDataType::STRING:
        LuaCheckStringValue(fname, L, 3);
        option.SetStringValue(lua_tostring(L, 3));
        opensn::Chi::log.Log() << "Solver:" << solver.TextName() << " option:" << option_name
                               << " set to " << option.StringValue() << ".";
        break;
      case VaryingDataType::BOOL:
        LuaCheckBoolValue(fname, L, 3);
        option.SetBoolValue(lua_toboolean(L, 3));
        opensn::Chi::log.Log() << "Solver:" << solver.TextName() << " option:" << option_name
                               << " set to " << ((option.BoolValue()) ? "true" : "false") << ".";
        break;
      case VaryingDataType::INTEGER:
        LuaCheckIntegerValue(fname, L, 3);
        option.SetIntegerValue(lua_tointeger(L, 3));
        opensn::Chi::log.Log() << "Solver:" << solver.TextName() << " option:" << option_name
                               << " set to " << option.IntegerValue() << ".";
        break;
      case VaryingDataType::FLOAT:
        LuaCheckNumberValue(fname, L, 3);
        option.SetFloatValue(lua_tonumber(L, 3));
        opensn::Chi::log.Log() << "Solver:" << solver.TextName() << " option:" << option_name
                               << " set to " << option.FloatValue() << ".";
        break;
    }
  }
  catch (const std::out_of_range& oor)
  {
    opensn::Chi::log.Log0Error() << fname << ": " << oor.what();
    throw oor;
  }

  return 0;
}

int
chiSolverGetName(lua_State* L)
{
  const std::string fname = "chiSolverGetName";
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  const auto& solver =
    opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  lua_pushstring(L, solver.TextName().c_str());

  return 1;
}

int
chiSolverGetFieldFunctionList(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError("chiGetFieldFunctionList", 1, num_args);

  // Getting solver
  const int solver_handle = lua_tonumber(L, 1);

  const auto& solver =
    opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  // Push up new table
  lua_newtable(L);
  for (size_t ff = 0; ff < solver.GetFieldFunctions().size(); ff++)
  {
    lua_pushinteger(L, static_cast<lua_Integer>(ff) + 1);
    int pff_count = -1;
    bool found = false;
    for (auto& pff : opensn::Chi::field_function_stack) // pff pointer to field func
    {
      ++pff_count;
      if (pff == solver.GetFieldFunctions()[ff])
      {
        lua_pushnumber(L, pff_count);
        found = true;
        break;
      }
    }

    if (not found)
      throw std::logic_error(fname + ": The solver specified has no "
                                     "field functions that match the global"
                                     " stack.");
    lua_settable(L, -3);
  }

  lua_pushinteger(L, static_cast<lua_Integer>(solver.GetFieldFunctions().size()));

  return 2;
}

int
chiSolverGetInfo(lua_State* L)
{
  const std::string fname = "chiSolverGetInfo";
  const int num_args = lua_gettop(L);

  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const size_t solver_handle = lua_tointeger(L, 1);

  const auto& solver =
    opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  ParameterBlock params;
  if (lua_isstring(L, 2)) params.AddParameter("name", std::string(lua_tostring(L, 2)));
  else if (lua_istable(L, 2))
    params = TableParserAsParameterBlock::ParseTable(L, 2);
  else
    ChiInvalidArgument("Argument 2 can only take a string or a table");

  const auto output_params = solver.GetInfo(params);

  PushParameterBlock(L, output_params);

  const int num_sub_params = static_cast<int>(output_params.NumParameters());

  return output_params.IsScalar() ? 1 : num_sub_params;
}

int
chiSolverSetProperties(lua_State* L)
{
  const std::string fname = "chiSolverSetProperties";
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const size_t solver_handle = lua_tointeger(L, 1);

  auto& solver = opensn::Chi::GetStackItem<Solver>(opensn::Chi::object_stack, solver_handle, fname);

  LuaCheckTableValue(fname, L, 2);
  auto property_block = TableParserAsParameterBlock::ParseTable(L, 2);

  solver.SetProperties(property_block);

  return 0;
}

} // namespace opensnlua
