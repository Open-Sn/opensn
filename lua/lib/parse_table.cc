// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/parse_table.h"
#include "framework/materials/material_property.h"
#include "framework/mesh/mesh_vector.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/parameters/parameter_block.h"
#include "framework/parameters/input_parameters.h"
#include "framework/math/functions/function.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/graphs/graph_partitioner.h"
#include "framework/post_processors/post_processor.h"
#include "framework/physics/solver.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/point_source/point_source.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/volumetric_source/volumetric_source.h"
#include "modules/linear_boltzmann_solvers/response_evaluator/response_evaluator.h"
#include "lua/lib/types.h"
#include "LuaBridge/LuaBridge.h"

namespace opensnlua
{

template <typename T>
T
CreateObject(lua_State* L)
{
  return luabridge::Stack<T>::get(L, lua_gettop(L)).value();
}

template <typename T>
std::shared_ptr<T>
CreateObjectPtr(lua_State* L)
{
  return luabridge::Stack<std::shared_ptr<T>>::get(L, lua_gettop(L)).value();
}

void
SetBlockParam(lua_State* L,
              const std::string& cls_name,
              const std::string& key,
              opensn::ParameterBlock& block)
{
  if (cls_name == "Vector3")
    block.AddParameter(key, CreateObject<opensn::Vector3>(L));
  else if (cls_name == "MeshContinuum")
    block.AddParameter(key, CreateObjectPtr<opensn::MeshContinuum>(L));
  else if (cls_name == "SurfaceMesh")
    block.AddParameter(key, CreateObjectPtr<opensn::SurfaceMesh>(L));
  // mesh generators
  else if (cls_name == "MeshGenerator")
    block.AddParameter(key, CreateObjectPtr<opensn::MeshGenerator>(L));
  else if (cls_name == "ExtruderMeshGenerator")
    block.AddParameter(key, CreateObjectPtr<opensn::MeshGenerator>(L));
  else if (cls_name == "OrthogonalMeshGenerator")
    block.AddParameter(key, CreateObjectPtr<opensn::MeshGenerator>(L));
  else if (cls_name == "FromFileMeshGenerator")
    block.AddParameter(key, CreateObjectPtr<opensn::MeshGenerator>(L));
  else if (cls_name == "SplitFileMeshGenerator")
    block.AddParameter(key, CreateObjectPtr<opensn::MeshGenerator>(L));
  // logical volumes
  else if (cls_name == "BooleanLogicalVolume")
    block.AddParameter(key, CreateObjectPtr<opensn::LogicalVolume>(L));
  else if (cls_name == "RCCLogicalVolume")
    block.AddParameter(key, CreateObjectPtr<opensn::LogicalVolume>(L));
  else if (cls_name == "RPPLogicalVolume")
    block.AddParameter(key, CreateObjectPtr<opensn::LogicalVolume>(L));
  else if (cls_name == "SphereLogicalVolume")
    block.AddParameter(key, CreateObjectPtr<opensn::LogicalVolume>(L));
  else if (cls_name == "SphereLogicalVolume")
    block.AddParameter(key, CreateObjectPtr<opensn::LogicalVolume>(L));
  else if (cls_name == "SurfaceMeshLogicalVolume")
    block.AddParameter(key, CreateObjectPtr<opensn::LogicalVolume>(L));
  //
  else if (cls_name == "KBAGraphPartitioner")
    block.AddParameter(key, CreateObjectPtr<opensn::GraphPartitioner>(L));
  else if (cls_name == "PETScGraphPartitioner")
    block.AddParameter(key, CreateObjectPtr<opensn::GraphPartitioner>(L));
  //
  else if (cls_name == "FieldFunctionGridBased")
    block.AddParameter(key, CreateObjectPtr<opensn::FieldFunction>(L));
  else if (cls_name == "FieldFunctionInterpolation")
    block.AddParameter(key, CreateObjectPtr<opensn::FieldFunctionInterpolation>(L));
  //
  else if (cls_name == "AngularQuadrature")
    block.AddParameter(key, CreateObjectPtr<opensn::AngularQuadrature>(L));
  else if (cls_name == "ProductQuadrature")
    block.AddParameter(key, CreateObjectPtr<opensn::AngularQuadrature>(L));
  else if (cls_name == "SLDFESQuadrature")
    block.AddParameter(key, CreateObjectPtr<opensn::AngularQuadrature>(L));
  //
  else if (cls_name == "MultiGroupXS")
    block.AddParameter(key, CreateObjectPtr<opensn::MaterialProperty>(L));
  else if (cls_name == "IsotropicMultiGroupSource")
    block.AddParameter(key, CreateObjectPtr<opensn::MaterialProperty>(L));
  //
  else if (cls_name == "AggregateNodalValuePostProcessor")
    block.AddParameter(key, CreateObjectPtr<opensn::PostProcessor>(L));
  //
  else if (cls_name == "PRKSolver")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "LBSSolver")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "DiscreteOrdinatesSolver")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "DiscreteOrdinatesCurvilinearSolver")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "NonLinearKEigen")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "PowerIterationKEigen")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "PowerIterationKEigenSCDSA")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "PowerIterationKEigenSMM")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "SteadyStateSolver")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  else if (cls_name == "DiffusionDFEMSolver")
    block.AddParameter(key, CreateObjectPtr<opensn::Solver>(L));
  //
  else if (cls_name == "PointSource")
    block.AddParameter(key, CreateObjectPtr<opensn::PointSource>(L));
  else if (cls_name == "VolumetricSource")
    block.AddParameter(key, CreateObjectPtr<opensn::VolumetricSource>(L));
  //
  else if (cls_name == "ScalarMaterialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "ScalarSpatialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "ScalarSpatialMaterialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "SpatialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "VectorSpatialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "VectorSpatialMaterialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "LuaScalarMaterialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "LuaScalarSpatialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "LuaScalarSpatialMaterialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "LuaSpatialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "LuaVectorSpatialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  else if (cls_name == "LuaVectorSpatialMaterialFunction")
    block.AddParameter(key, CreateObjectPtr<opensn::Function>(L));
  //
  else if (cls_name == "ResponseEvaluator")
    block.AddParameter(key, CreateObjectPtr<opensn::ResponseEvaluator>(L));
  else
    throw std::runtime_error("Unsupported type on Lua stack: " + cls_name);
}

void
ParseTableValues(lua_State* L, opensn::ParameterBlock& block, const std::string& key)
{
  auto v = luabridge::LuaRef::fromStack(L, -1);
  switch (v.type())
  {
    case LUA_TNIL:
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Encountered nil value assigned to key " + key);

    case LUA_TBOOLEAN:
      block.AddParameter(key, luabridge::Stack<bool>::get(L, -1).value());
      break;

    case LUA_TNUMBER:
      if (lua_isinteger(L, -1))
        block.AddParameter(key, luabridge::Stack<int64_t>::get(L, -1).value());
      else
        block.AddParameter(key, luabridge::Stack<double>::get(L, -1).value());
      break;

    case LUA_TSTRING:
      block.AddParameter(key, luabridge::Stack<std::string>::get(L, -1).value());
      break;

    case LUA_TTABLE:
    {
      auto v = luabridge::LuaRef::fromStack(L, lua_gettop(L));
      opensn::ParameterBlock new_block(key);
      ParseTableKeys(L, lua_gettop(L), new_block);
      block.AddParameter(new_block);
      break;
    }

    case LUA_TUSERDATA:
    {
      auto v = luabridge::LuaRef::fromStack(L, lua_gettop(L));
      auto cls_name = v.getClassName().value();
      SetBlockParam(L, cls_name, key, block);
      break;
    }

    default:
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Encountered unsupported value type " +
                             lua_typename(L, lua_type(L, -2)) + " for key " + key);
  }
}

void
ParseTableKeys(lua_State* L, int idx, opensn::ParameterBlock& block)
{
  bool number_key_encountered = false;
  bool string_key_encountered = false;
  int key_number_index = 0;

  lua_pushnil(L);
  while (lua_next(L, idx) != 0)
  {
    const int top = lua_gettop(L);
    if (lua_type(L, -2) == LUA_TSTRING)
    {
      if (number_key_encountered)
        throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                               ": Encountered mixed key types (string and number)");

      string_key_encountered = true;
      auto key = luabridge::Stack<std::string>::get(L, top - 1).value();
      ParseTableValues(L, block, key);
    } // if key is string

    if (lua_type(L, -2) == LUA_TNUMBER)
    {
      // If the key is a number then the following apply:
      // - This must be an array of items
      // - All the keys in the table must be numbers
      if (string_key_encountered)
        throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                               ": Encountered mixed key types (string and number)");

      if (block.GetType() != opensn::ParameterBlockType::ARRAY)
      {
        block.ChangeToArray();
      }

      number_key_encountered = true;
      const std::string key = std::to_string(key_number_index);
      ParseTableValues(L, block, key);
      ++key_number_index;
    }

    lua_pop(L, 1);
  }
}

luabridge::Result
PushParameterBlock(lua_State* L, const opensn::ParameterBlock& block, int level)
{
  switch (block.GetType())
  {
    case opensn::ParameterBlockType::BOOLEAN:
      return luabridge::push(L, block.GetValue<bool>());
    case opensn::ParameterBlockType::FLOAT:
      return luabridge::push(L, block.GetValue<double>());
    case opensn::ParameterBlockType::STRING:
      return luabridge::push(L, block.GetValue<std::string>());
    case opensn::ParameterBlockType::INTEGER:
      return luabridge::push(L, block.GetValue<int64_t>());
    case opensn::ParameterBlockType::ARRAY:
    {
      luabridge::Result result;
      if (level > 0)
        lua_newtable(L);
      const size_t num_params = block.GetNumParameters();
      for (size_t k = 0; k < num_params; ++k)
      {
        if (level > 0)
        {
          result = luabridge::push(L, k + 1);
          if (not result)
            return result;
        }
        result = PushParameterBlock(L, block.GetParam(k), level + 1);
        if (not result)
          return result;
        if (level > 0)
          lua_settable(L, -3);
      }
      return {};
    }
    case opensn::ParameterBlockType::BLOCK:
    {
      luabridge::Result result;
      if (level > 0)
        lua_newtable(L);
      const size_t num_params = block.GetNumParameters();
      for (size_t k = 0; k < num_params; ++k)
      {
        const auto& param = block.GetParam(k);
        if (level > 0)
        {
          result = luabridge::push(L, param.GetName());
          if (not result)
            return result;
        }
        result = PushParameterBlock(L, block.GetParam(k), level + 1);
        if (not result)
          return result;
        if (level > 0)
          lua_settable(L, -3);
      }
      break;
    }
    default:
      OpenSnLogicalError("Attempting to push unsupported ParameterBlockType to lua");
  }
  return {};
}

} // namespace opensnlua
