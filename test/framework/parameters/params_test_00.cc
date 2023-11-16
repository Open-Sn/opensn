#include "framework/parameters/parameter_block.h"

#include "framework/lua.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "framework/console/console.h"

using namespace opensn;

namespace chi_unit_tests
{
int chi_ParameterBlock_Test00(lua_State* L);

RegisterLuaFunction(chi_ParameterBlock_Test00, chi_unit_tests, chi_ParameterBlock_Test00);

int
chi_ParameterBlock_Test00(lua_State* L)
{
  opensn::Chi::log.Log() << "GOLD_BEGIN";
  const int num_args = lua_gettop(L);
  bool verbose = false;
  if (num_args >= 1) verbose = lua_toboolean(L, 1);

  if (verbose) opensn::Chi::log.Log() << "Hello world";

  if (num_args == 2)
  {
    if (lua_istable(L, 2))
    {
      opensn::Chi::log.Log() << "It is a block";
      const auto param_block = TableParserAsParameterBlock::ParseTable(L, 2);

      {
        std::string outstr;
        param_block.RecursiveDumpToString(outstr);

        opensn::Chi::log.Log() << outstr;
      }

      opensn::Chi::log.Log() << param_block.GetParamValue<std::string>("it_method");

      opensn::Chi::log.Log() << param_block.GetParam("sub1").GetParamValue<int>("ax_method");

      opensn::Chi::log.Log() << param_block.GetParamValue<double>("nl_abs_tol");

      opensn::Chi::log.Log() << (param_block.GetParamValue<bool>("enabled") ? "true" : "false");

      opensn::Chi::log.Log() << param_block.GetParamValue<size_t>("nl_max_its");

      opensn::Chi::log.Log() << "Has \"blocks\"?: " << param_block.GetParam("sub2").Has("blocks");

      opensn::Chi::log.Log() << "Num Parameters: "
                             << param_block.GetParam("sub2").GetParam("blocks").NumParameters();

      const auto vec = param_block.GetParam("sub2").GetParamVectorValue<int>("blocks");

      {
        std::stringstream outstr;
        for (auto val : vec)
          outstr << val << " ";
        opensn::Chi::log.Log() << outstr.str();
      }

      opensn::Chi::log.Log() << "Testing copy constructor";
      const auto& param_block2 = param_block;
      {
        std::string outstr;
        param_block2.RecursiveDumpToString(outstr);

        opensn::Chi::log.Log() << outstr;
      }

      opensn::Chi::log.Log() << "Testing move constructor";
      const ParameterBlock& param_block3(param_block2);
      {
        std::string outstr;
        param_block3.RecursiveDumpToString(outstr);

        opensn::Chi::log.Log() << outstr;
      }
    } // if table
  }   // if num_args == 2

  opensn::Chi::log.Log() << "Testing varying";
  {
    Varying v(12);
    opensn::Chi::log.Log() << "v(12)" << v.IntegerValue();
    v = true;
    opensn::Chi::log.Log() << "v(bool)" << v.BoolValue();
    opensn::Chi::log.Log() << "v(bool)" << v.GetValue<bool>();
  }
  {
    Varying v(12);
    opensn::Chi::log.Log() << "v(12)" << v.IntegerValue();
    v = 12.0;
    opensn::Chi::log.Log() << "v(12.0)" << v.FloatValue();
    opensn::Chi::log.Log() << "v(12.0)" << v.GetValue<double>();
    opensn::Chi::log.Log() << "v(12.0)" << v.GetValue<float>();
  }
  {
    Varying v(12.0);
    opensn::Chi::log.Log() << "v(12.0) bytesize" << v.ByteSize();
    opensn::Chi::log.Log() << "v(12.0)" << v.FloatValue();
    v = 12;
    opensn::Chi::log.Log() << "v(12)" << v.IntegerValue();
    opensn::Chi::log.Log() << "v(12)" << v.GetValue<int>();
    opensn::Chi::log.Log() << "v(12)" << v.GetValue<size_t>();
  }
  {
    Varying v(std::string("Hello"));
    opensn::Chi::log.Log() << "hello" << v.StringValue();
    opensn::Chi::log.Log() << "hello" << v.GetValue<std::string>();
  }
  opensn::Chi::log.Log() << "GOLD_END";
  return 0;
}

} // namespace chi_unit_tests
