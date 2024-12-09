#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/lib/console.h"

using namespace opensn;
using namespace opensnlua;

namespace unit_tests
{

void
ParameterBlock_Test00(const ParameterBlock& param_block)
{
  opensn::log.Log() << "GOLD_BEGIN";

  opensn::log.Log() << "Hello world";

  opensn::log.Log() << "It is a block";

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

  opensn::log.Log() << "Testing varying";
  {
    Varying v(12);
    opensn::log.Log() << "v(12)" << v.GetIntegerValue();
    v = true;
    opensn::log.Log() << "v(bool)" << v.GetBoolValue();
    opensn::log.Log() << "v(bool)" << v.GetValue<bool>();
  }
  {
    Varying v(12);
    opensn::log.Log() << "v(12)" << v.GetIntegerValue();
    v = 12.0;
    opensn::log.Log() << "v(12.0)" << v.GetFloatValue();
    opensn::log.Log() << "v(12.0)" << v.GetValue<double>();
    opensn::log.Log() << "v(12.0)" << v.GetValue<float>();
  }
  {
    Varying v(12.0);
    opensn::log.Log() << "v(12.0) bytesize" << v.GetByteSize();
    opensn::log.Log() << "v(12.0)" << v.GetFloatValue();
    v = 12;
    opensn::log.Log() << "v(12)" << v.GetIntegerValue();
    opensn::log.Log() << "v(12)" << v.GetValue<int>();
    opensn::log.Log() << "v(12)" << v.GetValue<size_t>();
  }
  {
    Varying v(std::string("Hello"));
    opensn::log.Log() << "hello" << v.GetStringValue();
    opensn::log.Log() << "hello" << v.GetValue<std::string>();
  }
  opensn::log.Log() << "GOLD_END";
}

BIND_FUNCTION(unit_tests, ParameterBlock_Test00);

} //  namespace unit_tests
