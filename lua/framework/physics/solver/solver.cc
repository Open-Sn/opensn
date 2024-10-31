// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/physics/solver/solver.h"
#include "lua/framework/console/console.h"
#include "framework/physics/solver.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/event_system/physics_event_publisher.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(SolverCreate, solver, Create);
RegisterLuaFunctionInNamespace(SolverInitialize, solver, Initialize);
RegisterLuaFunctionInNamespace(SolverExecute, solver, Execute);
RegisterLuaFunctionInNamespace(SolverStep, solver, Step);
RegisterLuaFunctionInNamespace(SolverAdvance, solver, Advance);
RegisterLuaFunctionInNamespace(SolverSetBasicOption, solver, SetBasicOption);
RegisterLuaFunctionInNamespace(SolverGetName, solver, GetName);
RegisterLuaFunctionInNamespace(SolverGetFieldFunctionList, solver, GetFieldFunctionList);
RegisterLuaFunctionInNamespace(SolverGetInfo, solver, GetInfo);
RegisterLuaFunctionInNamespace(SolverSetProperties, solver, SetProperties);

int
SolverCreate(lua_State* L)
{
  const std::string fname = "solver.Create";
  LuaCheckArgs<ParameterBlock>(L, fname);

  const auto params = LuaArg<ParameterBlock>(L, 1);
  const auto& object_maker = ObjectFactory::GetInstance();
  const size_t handle = object_maker.MakeRegisteredObject(params);
  return LuaReturn(L, handle);
}

int
SolverInitialize(lua_State* L)
{
  const std::string fname = "solver.Initialize";
  LuaCheckArgs<size_t>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);

  auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);

  PhysicsEventPublisher::GetInstance().SolverInitialize(solver);

  return LuaReturn(L);
}

int
SolverExecute(lua_State* L)
{
  const std::string fname = "solver.Execute";
  LuaCheckArgs<size_t>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);

  auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);

  PhysicsEventPublisher::GetInstance().SolverExecute(solver);

  return LuaReturn(L);
}

int
SolverStep(lua_State* L)
{
  const std::string fname = "solver.Step";
  LuaCheckArgs<size_t>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);

  auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);

  PhysicsEventPublisher::GetInstance().SolverStep(solver);

  return LuaReturn(L);
}

int
SolverAdvance(lua_State* L)
{
  const std::string fname = "solver.Advance";
  LuaCheckArgs<size_t>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);

  auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);

  PhysicsEventPublisher::GetInstance().SolverAdvance(solver);

  return LuaReturn(L);
}

int
SolverSetBasicOption(lua_State* L)
{
  const std::string fname = "solver.SetBasicOption";
  LuaCheckArgs<int, std::string>(L, fname);

  const auto solver_handle = LuaArg<int>(L, 1);
  const auto option_name = LuaArg<std::string>(L, 2);
  auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);

  try
  {
    auto& option = solver.GetBasicOptions()[option_name];

    switch (option.Type())
    {
      case VaryingDataType::VOID:
      case VaryingDataType::ARBITRARY_BYTES:
        throw std::logic_error("Solver:" + solver.Name() + " option:" + option_name +
                               " is of invalid type."
                               " This indicates an implementation problem.");
      case VaryingDataType::STRING:
        LuaCheckArgs<int, std::string, std::string>(L, fname);
        option.SetStringValue(LuaArg<std::string>(L, 3));
        opensn::log.Log() << "Solver:" << solver.Name() << " option:" << option_name << " set to "
                          << option.StringValue() << ".";
        break;
      case VaryingDataType::BOOL:
        LuaCheckArgs<int, std::string, bool>(L, fname);
        option.SetBoolValue(LuaArg<bool>(L, 3));
        opensn::log.Log() << "Solver:" << solver.Name() << " option:" << option_name << " set to "
                          << ((option.BoolValue()) ? "true" : "false") << ".";
        break;
      case VaryingDataType::INTEGER:
        LuaCheckArgs<int, std::string, int>(L, fname);
        option.SetIntegerValue(LuaArg<int>(L, 3));
        opensn::log.Log() << "Solver:" << solver.Name() << " option:" << option_name << " set to "
                          << option.IntegerValue() << ".";
        break;
      case VaryingDataType::FLOAT:
        LuaCheckArgs<int, std::string, double>(L, fname);
        option.SetFloatValue(LuaArg<double>(L, 3));
        opensn::log.Log() << "Solver:" << solver.Name() << " option:" << option_name << " set to "
                          << option.FloatValue() << ".";
        break;
    }
  }
  catch (const std::out_of_range& oor)
  {
    opensn::log.Log0Error() << fname << ": " << oor.what();
    throw oor;
  }

  return LuaReturn(L);
}

int
SolverGetName(lua_State* L)
{
  const std::string fname = "solver.GetName";
  LuaCheckArgs<int>(L, fname);

  const auto solver_handle = LuaArg<int>(L, 1);
  const auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);
  auto name = solver.Name();
  return LuaReturn(L, name);
}

int
SolverGetFieldFunctionList(lua_State* L)
{
  const std::string fname = "solver.GetFieldFunctionList";
  LuaCheckArgs<size_t>(L, fname);

  // Getting solver
  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);

  std::vector<size_t> ff_handles;
  for (size_t ff = 0; ff < solver.GetFieldFunctions().size(); ++ff)
  {
    int pff_count = -1;
    bool found = false;
    for (auto& pff : opensn::field_function_stack) // pff pointer to field func
    {
      ++pff_count;
      if (pff == solver.GetFieldFunctions()[ff])
      {
        ff_handles.push_back(pff_count);
        found = true;
        break;
      }
    }
    if (not found)
      throw std::logic_error(fname + ": The solver specified has no "
                                     "field functions that match the global"
                                     " stack.");
  }
  return LuaReturn(L, ff_handles, ff_handles.size());
}

int
SolverGetInfo(lua_State* L)
{
  const std::string fname = "solver.GetInfo";
  LuaCheckArgs<size_t>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);

  ParameterBlock params;
  if (lua_isstring(L, 2))
  {
    LuaCheckArgs<size_t, std::string>(L, fname);
    params.AddParameter("name", LuaArg<std::string>(L, 2));
  }
  else if (lua_istable(L, 2))
    params = LuaArg<ParameterBlock>(L, 2);
  else
    OpenSnInvalidArgument("Argument 2 can only take a string or a table");

  const auto output_params = solver.GetInfo(params);

  LuaPush(L, output_params);

  const int num_sub_params = static_cast<int>(output_params.NumParameters());

  return output_params.IsScalar() ? 1 : num_sub_params;
}

int
SolverSetProperties(lua_State* L)
{
  const std::string fname = "solver.SetProperties";
  LuaCheckArgs<size_t, ParameterBlock>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);
  auto& solver = opensn::GetStackItem<Solver>(opensn::object_stack, solver_handle, fname);
  auto property_block = LuaArg<ParameterBlock>(L, 2);

  solver.SetProperties(property_block);

  return LuaReturn(L);
}

} // namespace opensnlua
