#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "framework/console/chi_console.h"

#include "framework/utils/chi_timer.h"

namespace chi_unit_tests
{

chi::ParameterBlock LogTimingInfoTest(const chi::InputParameters&);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/LogTimingInfoTest,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/LogTimingInfoTest);

chi::ParameterBlock
LogTimingInfoTest(const chi::InputParameters&)
{
  Chi::log.Log() << "LogTiming test";
  auto& t_main = Chi::log.CreateTimingBlock("Timing_Main");
  t_main.TimeSectionBegin();
  {
    // Some overhead
    chi::Sleep(std::chrono::milliseconds(200));

    auto& t_1 = Chi::log.CreateOrGetTimingBlock("Part1", "Timing_Main");
    t_1.TimeSectionBegin();
    {
      chi::Sleep(std::chrono::milliseconds(300));
    }
    t_1.TimeSectionEnd();

    auto& t_2 = Chi::log.CreateOrGetTimingBlock("Part2", "Timing_Main");
    t_2.TimeSectionBegin();
    {
      chi::Sleep(std::chrono::milliseconds(333));
    }
    t_2.TimeSectionEnd();

    auto& t_3 = Chi::log.GetTimingBlock("Part2");
    t_3.TimeSectionBegin();
    {
      chi::Sleep(std::chrono::milliseconds(123));
    }
    t_3.TimeSectionEnd();
  }
  t_main.TimeSectionEnd();

  Chi::log.Log() << Chi::log.GetTimingBlock("ChiTech").MakeGraphString();

  return chi::ParameterBlock{};
}

} // namespace chi_unit_tests
