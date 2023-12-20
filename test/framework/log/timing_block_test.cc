#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "framework/console/console.h"

#include "framework/utils/timer.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock LogTimingInfoTest(const InputParameters&);

RegisterWrapperFunction(chi_unit_tests, LogTimingInfoTest, nullptr, LogTimingInfoTest);

ParameterBlock
LogTimingInfoTest(const InputParameters&)
{
  opensn::Chi::log.Log() << "LogTiming test";
  auto& t_main = opensn::Chi::log.CreateTimingBlock("Timing_Main");
  t_main.TimeSectionBegin();
  {
    // Some overhead
    opensn::Sleep(std::chrono::milliseconds(200));

    auto& t_1 = opensn::Chi::log.CreateOrGetTimingBlock("Part1", "Timing_Main");
    t_1.TimeSectionBegin();
    {
      opensn::Sleep(std::chrono::milliseconds(300));
    }
    t_1.TimeSectionEnd();

    auto& t_2 = opensn::Chi::log.CreateOrGetTimingBlock("Part2", "Timing_Main");
    t_2.TimeSectionBegin();
    {
      opensn::Sleep(std::chrono::milliseconds(333));
    }
    t_2.TimeSectionEnd();

    auto& t_3 = opensn::Chi::log.GetTimingBlock("Part2");
    t_3.TimeSectionBegin();
    {
      opensn::Sleep(std::chrono::milliseconds(123));
    }
    t_3.TimeSectionEnd();
  }
  t_main.TimeSectionEnd();

  opensn::Chi::log.Log() << opensn::Chi::log.GetTimingBlock("ChiTech").MakeGraphString();

  return ParameterBlock{};
}

} //  namespace unit_tests
