#include "framework/materials/multi_group_xs/xsfile.h"
#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

using namespace opensn;

namespace
{

std::string
MakeTempXSFile()
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  auto path =
    std::filesystem::temp_directory_path() / ("opensn_xsfile_test_" + std::to_string(now) + ".xs");

  std::ofstream out(path);
  out << "NUM_GROUPS 2\n";
  out << "GROUP_STRUCTURE_BEGIN\n";
  out << "1.0\n";
  out << "0.1\n";
  out << "0.01\n";
  out << "GROUP_STRUCTURE_END\n";
  out << "SIGMA_T_BEGIN\n";
  out << "0 1.0\n";
  out << "1 2.0\n";
  out << "SIGMA_T_END\n";
  out << "SIGMA_A_BEGIN\n";
  out << "0 0.5\n";
  out << "1 1.5\n";
  out << "SIGMA_A_END\n";
  out.close();

  return path.string();
}

} // namespace

TEST(XSFileTest, ReadMinimalFile)
{
  const std::string fname = MakeTempXSFile();
  XSFile xs(fname);
  xs.Read();

  EXPECT_EQ(xs.num_groups_, 2u);
  ASSERT_EQ(xs.e_bounds_.size(), 3u);
  EXPECT_NEAR(xs.e_bounds_[0], 1.0, 1e-12);
  EXPECT_NEAR(xs.e_bounds_[1], 0.1, 1e-12);
  EXPECT_NEAR(xs.e_bounds_[2], 0.01, 1e-12);

  ASSERT_EQ(xs.sigma_t_.size(), 2u);
  EXPECT_NEAR(xs.sigma_t_[0], 1.0, 1e-12);
  EXPECT_NEAR(xs.sigma_t_[1], 2.0, 1e-12);

  ASSERT_EQ(xs.sigma_a_.size(), 2u);
  EXPECT_NEAR(xs.sigma_a_[0], 0.5, 1e-12);
  EXPECT_NEAR(xs.sigma_a_[1], 1.5, 1e-12);

  std::filesystem::remove(fname);
}
