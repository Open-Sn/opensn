// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/functions.h"
#include "framework/runtime.h"
#include "lua/lib/console.h"
#include "lua/lib/types.h"
#include "LuaBridge/LuaBridge.h"
#include <lua.h>

using namespace opensn;

namespace opensnlua
{

OpenSnRegisterObjectInNamespace(opensn, LuaScalarMaterialFunction);

InputParameters
LuaScalarMaterialFunction::GetInputParameters()
{
  InputParameters params = ScalarMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaScalarMaterialFunction>
LuaScalarMaterialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaScalarMaterialFunction>("opensn::LuaScalarMaterialFunction", params);
}

LuaScalarMaterialFunction::LuaScalarMaterialFunction(const InputParameters& params)
  : ScalarMaterialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

double
LuaScalarMaterialFunction::Evaluate(double val, int mat_id) const
{
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  return fn(val, mat_id)[0];
}

//

OpenSnRegisterObjectInNamespace(opensn, LuaScalarSpatialFunction);

InputParameters
LuaScalarSpatialFunction::GetInputParameters()
{
  InputParameters params = ScalarSpatialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaScalarSpatialFunction>
LuaScalarSpatialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaScalarSpatialFunction>("opensn::LuaScalarSpatialFunction", params);
}

LuaScalarSpatialFunction::LuaScalarSpatialFunction(const InputParameters& params)
  : ScalarSpatialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

double
LuaScalarSpatialFunction::Evaluate(const opensn::Vector3& xyz) const
{
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  return fn(xyz)[0];
}

//

OpenSnRegisterObjectInNamespace(opensn, LuaScalarSpatialMaterialFunction);

InputParameters
LuaScalarSpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = ScalarSpatialMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaScalarSpatialMaterialFunction>
LuaScalarSpatialMaterialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaScalarSpatialMaterialFunction>(
    "opensn::LuaScalarSpatialMaterialFunction", params);
}

LuaScalarSpatialMaterialFunction::LuaScalarSpatialMaterialFunction(const InputParameters& params)
  : ScalarSpatialMaterialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

double
LuaScalarSpatialMaterialFunction::Evaluate(int mat_id, const opensn::Vector3& xyz) const
{
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  auto res = fn(mat_id, xyz);
  return res[0];
}

//

OpenSnRegisterObjectInNamespace(opensn, LuaVectorSpatialFunction);

InputParameters
LuaVectorSpatialFunction::GetInputParameters()
{
  InputParameters params = VectorSpatialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaVectorSpatialFunction>
LuaVectorSpatialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaVectorSpatialFunction>("opensn::LuaVectorSpatialFunction", params);
}

LuaVectorSpatialFunction::LuaVectorSpatialFunction(const opensn::InputParameters& params)
  : opensn::VectorSpatialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

std::vector<double>
LuaVectorSpatialFunction::Evaluate(const opensn::Vector3& xyz, int num_components) const
{
  // Load lua function
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  auto lua_return = fn(xyz);
  OpenSnLogicalErrorIf(lua_return.hasFailed(), "Call to " + function_name_ + " failed.");
  auto tab = lua_return[0];
  OpenSnLogicalErrorIf(not tab.isTable(),
                       "The return value from " + function_name_ + " must be a table.");
  // Check return value
  OpenSnLogicalErrorIf(tab.length() != num_components,
                       "Call to lua function " + function_name_ + " returned a vector of size " +
                         std::to_string(tab.length()) +
                         ", which is not the same as the number of groups " +
                         std::to_string(num_components) + ".");
  std::vector<double> ret(tab.length());
  for (auto i = 0; i < tab.length(); i++)
    ret[i] = tab[i + 1];
  return ret;
}

//

OpenSnRegisterObjectInNamespace(opensn, LuaVectorSpatialMaterialFunction);

InputParameters
LuaVectorSpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = VectorSpatialMaterialFunction::GetInputParameters();
  params.AddRequiredParameter<std::string>("function_name", "Name of the lua function");
  return params;
}

std::shared_ptr<LuaVectorSpatialMaterialFunction>
LuaVectorSpatialMaterialFunction::Create(const opensn::ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<LuaVectorSpatialMaterialFunction>(
    "opensn::LuaVectorSpatialMaterialFunction", params);
}

LuaVectorSpatialMaterialFunction::LuaVectorSpatialMaterialFunction(const InputParameters& params)
  : opensn::VectorSpatialMaterialFunction(params),
    function_name_(params.GetParamValue<std::string>("function_name"))
{
}

std::vector<double>
LuaVectorSpatialMaterialFunction::Evaluate(const opensn::Vector3& xyz,
                                           int mat_id,
                                           int num_components) const
{
  // Load lua function
  lua_State* L = console.GetConsoleState();
  auto fn = luabridge::getGlobal(L, function_name_.c_str());
  auto lua_return = fn(xyz, mat_id);

  // Check return value
  OpenSnLogicalErrorIf(lua_return.size() != num_components,
                       "Call to lua function " + function_name_ + " returned a vector of size " +
                         std::to_string(lua_return.size()) +
                         ", which is not the same as the number of groups " +
                         std::to_string(num_components) + ".");
  std::vector<double> ret(lua_return.size());
  for (auto i = 0; i < lua_return.size(); i++)
    ret[i] = lua_return[i];
  return ret;
}

} // namespace opensnlua
