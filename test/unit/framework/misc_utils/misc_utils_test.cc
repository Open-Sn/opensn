#include "framework/utils/utils.h"
#include <gmock/gmock.h>

using namespace opensn;

TEST(MiscUtilsTest, Progress)
{
  const unsigned int I = 4;
  const size_t N = 39;

  std::stringstream progress;
  for (size_t i = 0; i < N; ++i)
  {
    progress << PrintIterationProgress(i, N, I);
  }

  EXPECT_EQ(progress.str(), "  25.00  50.00  75.00 100.00");
}
