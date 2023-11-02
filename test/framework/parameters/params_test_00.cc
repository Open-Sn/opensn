#include"parameters/parameter_block.h"

#include "chi_lua.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"

namespace chi_unit_tests
{
int chi_ParameterBlock_Test00(lua_State* L);

RegisterLuaFunction(/*function_ptr=*/chi_ParameterBlock_Test00,
                    /*namespace_name=*/chi_unit_tests,
                    /*func_name=*/chi_ParameterBlock_Test00);

int chi_ParameterBlock_Test00(lua_State* L)
{
  Chi::log.Log() << "GOLD_BEGIN";
  const int num_args = lua_gettop(L);
  bool verbose = false;
  if (num_args >= 1)
    verbose = lua_toboolean(L,1);

  if (verbose) Chi::log.Log() << "Hello world";

  if (num_args == 2)
  {
    if (lua_istable(L, 2))
    {
      Chi::log.Log() << "It is a block";
      const auto param_block =
        chi_lua::TableParserAsParameterBlock::ParseTable(L, 2);

      {
        std::string outstr;
        param_block.RecursiveDumpToString(outstr);

        Chi::log.Log() << outstr;
      }

      Chi::log.Log() << param_block.GetParamValue<std::string>("it_method");

      Chi::log.Log() << param_block.GetParam("sub1").GetParamValue<int>(
        "ax_method");

      Chi::log.Log() << param_block.GetParamValue<double>("nl_abs_tol");

      Chi::log.Log() << (param_block.GetParamValue<bool>("enabled") ? "true"
                                                                    : "false");

      Chi::log.Log() << param_block.GetParamValue<size_t>("nl_max_its");

      Chi::log.Log() << "Has \"blocks\"?: "
                     << param_block.GetParam("sub2").Has("blocks");

      Chi::log.Log()
        << "Num Parameters: "
        << param_block.GetParam("sub2").GetParam("blocks").NumParameters();

      const auto vec =
        param_block.GetParam("sub2").GetParamVectorValue<int>("blocks");

      {
        std::stringstream outstr;
        for (auto val : vec)
          outstr << val << " ";
        Chi::log.Log() << outstr.str();
      }

      Chi::log.Log() << "Testing copy constructor";
      const auto& param_block2 = param_block;
      {
        std::string outstr;
        param_block2.RecursiveDumpToString(outstr);

        Chi::log.Log() << outstr;
      }

      Chi::log.Log() << "Testing move constructor";
      const chi::ParameterBlock& param_block3(param_block2);
      {
        std::string outstr;
        param_block3.RecursiveDumpToString(outstr);

        Chi::log.Log() << outstr;
      }
    }//if table
  }//if num_args == 2

  Chi::log.Log() << "Testing varying";
  {
    chi_data_types::Varying v(12);
    Chi::log.Log() << "v(12)" << v.IntegerValue();
    v = true;
    Chi::log.Log() << "v(bool)" << v.BoolValue();
    Chi::log.Log() << "v(bool)" << v.GetValue<bool>();
  }
  {
    chi_data_types::Varying v(12);
    Chi::log.Log() << "v(12)" << v.IntegerValue();
    v = 12.0;
    Chi::log.Log() << "v(12.0)" << v.FloatValue();
    Chi::log.Log() << "v(12.0)" << v.GetValue<double>();
    Chi::log.Log() << "v(12.0)" << v.GetValue<float>();
  }
  {
    chi_data_types::Varying v(12.0);
    Chi::log.Log() << "v(12.0) bytesize" << v.ByteSize();
    Chi::log.Log() << "v(12.0)" << v.FloatValue();
    v = 12;
    Chi::log.Log() << "v(12)" << v.IntegerValue();
    Chi::log.Log() << "v(12)" << v.GetValue<int>();
    Chi::log.Log() << "v(12)" << v.GetValue<size_t>();
  }
  {
    chi_data_types::Varying v(std::string("Hello"));
    Chi::log.Log() << "hello" << v.StringValue();
    Chi::log.Log() << "hello" << v.GetValue<std::string>();
  }
  Chi::log.Log() << "GOLD_END";
  return 0;
}

}//namespace chi_unit_tests
