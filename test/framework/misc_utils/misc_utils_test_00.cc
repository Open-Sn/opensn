#include "framework/utils/utils.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "lua/lib/console.h"

using namespace opensn;

namespace unit_tests
{

void
misc_utils_Test00()
{
  opensn::log.Log() << "GOLD_BEGIN";
  opensn::log.Log() << "Testing misc_utils::PrintIterationProgress\n";

  const unsigned int I = 4;
  const size_t N = 39;

  std::stringstream progress;
  for (size_t i = 0; i < N; ++i)
  {
    progress << PrintIterationProgress(i, N, I);
  }

  opensn::log.Log() << progress.str();

  opensn::log.Log() << "GOLD_END";
}

BIND_FUNCTION(unit_tests, misc_utils_Test00);

} //  namespace unit_tests
