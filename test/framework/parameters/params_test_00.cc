#include "framework/parameters/parameter_block.h"

#include "lua/framework/lua.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "lua/framework/console/console.h"

using namespace opensn;
using namespace opensnlua;

namespace unit_tests
{
int ParameterBlock_Test00(lua_State* L);

RegisterLuaFunctionNamespace(ParameterBlock_Test00, unit_tests, ParameterBlock_Test00);

int
ParameterBlock_Test00(lua_State* L)
{
  opensn::log.Log() << "GOLD_BEGIN";
  const int num_args = lua_gettop(L);
  auto verbose = LuaArgOptional<bool>(L, 1, false);

  if (verbose)
    opensn::log.Log() << "Hello world";

  if (num_args == 2)
  {
    if (lua_istable(L, 2))
    {
      opensn::log.Log() << "It is a block";
      const auto param_block = opensnlua::TableParserAsParameterBlock::ParseTable(L, 2);

      {
        std::string outstr;
        param_block.RecursiveDumpToString(outstr);

        opensn::log.Log() << outstr;
      }

      opensn::log.Log() << param_block.GetParamValue<std::string>("it_method");

      opensn::log.Log() << param_block.GetParam("sub1").GetParamValue<int>("ax_method");

      opensn::log.Log() << param_block.GetParamValue<double>("nl_abs_tol");

      opensn::log.Log() << (param_block.GetParamValue<bool>("enabled") ? "true" : "false");

      opensn::log.Log() << param_block.GetParamValue<size_t>("nl_max_its");

      opensn::log.Log() << "Has \"blocks\"?: " << param_block.GetParam("sub2").Has("blocks");

      opensn::log.Log() << "Num Parameters: "
                        << param_block.GetParam("sub2").GetParam("blocks").NumParameters();

      const auto vec = param_block.GetParam("sub2").GetParamVectorValue<int>("blocks");

      {
        std::stringstream outstr;
        for (auto val : vec)
          outstr << val << " ";
        opensn::log.Log() << outstr.str();
      }

      opensn::log.Log() << "Testing copy constructor";
      const auto& param_block2 = param_block;
      {
        std::string outstr;
        param_block2.RecursiveDumpToString(outstr);

        opensn::log.Log() << outstr;
      }

      opensn::log.Log() << "Testing move constructor";
      const ParameterBlock& param_block3(param_block2);
      {
        std::string outstr;
        param_block3.RecursiveDumpToString(outstr);

        opensn::log.Log() << outstr;
      }
    } // if table
  }   // if num_args == 2

  opensn::log.Log() << "Testing varying";
  {
    Varying v(12);
    opensn::log.Log() << "v(12)" << v.IntegerValue();
    v = true;
    opensn::log.Log() << "v(bool)" << v.BoolValue();
    opensn::log.Log() << "v(bool)" << v.GetValue<bool>();
  }
  {
    Varying v(12);
    opensn::log.Log() << "v(12)" << v.IntegerValue();
    v = 12.0;
    opensn::log.Log() << "v(12.0)" << v.FloatValue();
    opensn::log.Log() << "v(12.0)" << v.GetValue<double>();
    opensn::log.Log() << "v(12.0)" << v.GetValue<float>();
  }
  {
    Varying v(12.0);
    opensn::log.Log() << "v(12.0) bytesize" << v.ByteSize();
    opensn::log.Log() << "v(12.0)" << v.FloatValue();
    v = 12;
    opensn::log.Log() << "v(12)" << v.IntegerValue();
    opensn::log.Log() << "v(12)" << v.GetValue<int>();
    opensn::log.Log() << "v(12)" << v.GetValue<size_t>();
  }
  {
    Varying v(std::string("Hello"));
    opensn::log.Log() << "hello" << v.StringValue();
    opensn::log.Log() << "hello" << v.GetValue<std::string>();
  }
  opensn::log.Log() << "GOLD_END";
  return 0;
}

} //  namespace unit_tests
