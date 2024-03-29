#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "lua/framework/console/console.h"

#include "framework/utils/timer.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock LogTimingInfoTest(const InputParameters&);

RegisterWrapperFunctionNamespace(unit_tests, LogTimingInfoTest, nullptr, LogTimingInfoTest);

ParameterBlock
LogTimingInfoTest(const InputParameters&)
{
  opensn::log.Log() << "LogTiming test";
  auto& t_main = opensn::log.CreateTimingBlock("Timing_Main");
  t_main.TimeSectionBegin();
  {
    // Some overhead
    opensn::Sleep(std::chrono::milliseconds(200));

    auto& t_1 = opensn::log.CreateOrGetTimingBlock("Part1", "Timing_Main");
    t_1.TimeSectionBegin();
    {
      opensn::Sleep(std::chrono::milliseconds(300));
    }
    t_1.TimeSectionEnd();

    auto& t_2 = opensn::log.CreateOrGetTimingBlock("Part2", "Timing_Main");
    t_2.TimeSectionBegin();
    {
      opensn::Sleep(std::chrono::milliseconds(333));
    }
    t_2.TimeSectionEnd();

    auto& t_3 = opensn::log.GetTimingBlock("Part2");
    t_3.TimeSectionBegin();
    {
      opensn::Sleep(std::chrono::milliseconds(123));
    }
    t_3.TimeSectionEnd();
  }
  t_main.TimeSectionEnd();

  opensn::log.Log() << opensn::log.GetTimingBlock(opensn::name).MakeGraphString();

  return ParameterBlock{};
}

} //  namespace unit_tests
