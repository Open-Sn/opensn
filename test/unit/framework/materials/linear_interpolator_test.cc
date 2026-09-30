#include "framework/materials/interpolator/cartesian_grid.h"
#include "framework/materials/interpolator/interpolator.h"
#include "framework/materials/interpolator/linear_interpolator.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"

#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <memory>
#include <vector>

using namespace opensn;

namespace
{

double
OneGroupTotal(double x, double y, double z)
{
  return 2.0 + 0.2 * x - 0.1 * y + 0.05 * z;
}

double
OneGroupScatter(double x, double y, double z)
{
  return 0.4 + 0.03 * x + 0.02 * y - 0.01 * z;
}

std::shared_ptr<MultiGroupXS>
MakeSimpleOneGroupXS(double x, double y, double z)
{
  const double sigma_t = OneGroupTotal(x, y, z);
  const double sigma_s = OneGroupScatter(x, y, z);
  return std::make_shared<MultiGroupXS>(
    MultiGroupXS::CreateSimpleOneGroup(sigma_t, sigma_s / sigma_t));
}

double
TwoGroupTotal(std::size_t g, double x, double y, double z)
{
  return 1.5 + 0.7 * static_cast<double>(g + 1) + 0.2 * x - 0.15 * y + 0.08 * z;
}

double
TwoGroupAbsorption(std::size_t g, double x, double y, double z)
{
  return 0.4 + 0.2 * static_cast<double>(g + 1) + 0.05 * x + 0.03 * y - 0.02 * z;
}

double
TwoGroupScatter(std::size_t ell, std::size_t row, std::size_t col, double x, double y, double z)
{
  return 0.1 * static_cast<double>(ell + 1) + 0.2 * static_cast<double>(row + 1) +
         0.3 * static_cast<double>(col + 1) + 0.04 * x - 0.01 * y + 0.02 * z;
}

std::shared_ptr<MultiGroupXS>
MakeTwoGroupP2XS(double x, double y, double z)
{
  auto xs = std::make_shared<MultiGroupXS>(2, 2);

  xs->GetSigmaTotal() = {TwoGroupTotal(0, x, y, z), TwoGroupTotal(1, x, y, z)};
  xs->GetSigmaAbsorption() = {TwoGroupAbsorption(0, x, y, z), TwoGroupAbsorption(1, x, y, z)};

  auto& transfer_matrices = xs->GetTransferMatrices();
  transfer_matrices.reserve(3);
  for (std::size_t ell = 0; ell < 3; ++ell)
  {
    SparseMatrix mat(2, 2);
    for (std::size_t row = 0; row < 2; ++row)
      for (std::size_t col = 0; col < 2; ++col)
        mat.Insert(row, col, TwoGroupScatter(ell, row, col, x, y, z));

    transfer_matrices.push_back(mat);
  }

  return xs;
}

} // namespace

TEST(LinearInterpolatorTest, EvaluateInterpolatesOneGroupTotalAbsorptionAndScatter)
{
  const std::vector<std::vector<double>> grid_data = {{0.0, 1.5}, {-1.0, 2.0}, {10.0, 13.0}};
  const CartesianGrid grid(grid_data);

  std::vector<std::shared_ptr<MultiGroupXS>> xs;
  xs.reserve(8);
  for (const double x : grid_data[0])
    for (const double y : grid_data[1])
      for (const double z : grid_data[2])
        xs.push_back(MakeSimpleOneGroupXS(x, y, z));

  LinearInterpolator interpolator(grid, xs, XSType::Total | XSType::Absorption | XSType::Transfer);

  std::vector<double> state_point = {0.63, 0.4, 11.7};
  const auto result = interpolator.Evaluate(state_point);

  const double expected_total = OneGroupTotal(state_point[0], state_point[1], state_point[2]);
  const double expected_scatter = OneGroupScatter(state_point[0], state_point[1], state_point[2]);
  const double expected_absorption = expected_total - expected_scatter;

  ASSERT_EQ(result.GetNumGroups(), 1u);
  ASSERT_EQ(result.GetSigmaTotal().size(), 1u);
  ASSERT_EQ(result.GetSigmaAbsorption().size(), 1u);
  ASSERT_EQ(result.GetTransferMatrices().size(), 1u);

  EXPECT_NEAR(result.GetSigmaTotal()[0], expected_total, 1.0e-12);
  EXPECT_NEAR(result.GetSigmaAbsorption()[0], expected_absorption, 1.0e-12);
  EXPECT_NEAR(result.GetTransferMatrix(0).GetValueIJ(0, 0), expected_scatter, 1.0e-12);
}

TEST(LinearInterpolatorTest, EvaluateInterpolatesTwoGroupP2TotalAbsorptionAndScatter)
{
  const std::vector<std::vector<double>> grid_data = {{0.0, 2.0}, {1.0, 3.0}, {-2.0, 1.0}};
  const CartesianGrid grid(grid_data);

  std::vector<std::shared_ptr<MultiGroupXS>> xs;
  xs.reserve(8);
  for (const double x : grid_data[0])
    for (const double y : grid_data[1])
      for (const double z : grid_data[2])
        xs.push_back(MakeTwoGroupP2XS(x, y, z));

  LinearInterpolator interpolator(grid, xs, XSType::Total | XSType::Absorption | XSType::Transfer);

  std::vector<double> state_point = {0.7, 2.2, -0.4};
  const auto result = interpolator.Evaluate(state_point);

  ASSERT_EQ(result.GetNumGroups(), 2u);
  ASSERT_EQ(result.GetScatteringOrder(), 2u);
  ASSERT_EQ(result.GetSigmaTotal().size(), 2u);
  ASSERT_EQ(result.GetSigmaAbsorption().size(), 2u);
  ASSERT_EQ(result.GetTransferMatrices().size(), 3u);

  for (std::size_t g = 0; g < 2; ++g)
  {
    EXPECT_NEAR(result.GetSigmaTotal()[g],
                TwoGroupTotal(g, state_point[0], state_point[1], state_point[2]),
                1.0e-12);
    EXPECT_NEAR(result.GetSigmaAbsorption()[g],
                TwoGroupAbsorption(g, state_point[0], state_point[1], state_point[2]),
                1.0e-12);
  }

  for (std::size_t ell = 0; ell < 3; ++ell)
    for (std::size_t row = 0; row < 2; ++row)
      for (std::size_t col = 0; col < 2; ++col)
        EXPECT_NEAR(result.GetTransferMatrix(ell).GetValueIJ(row, col),
                    TwoGroupScatter(ell, row, col, state_point[0], state_point[1], state_point[2]),
                    1.0e-12);
}

TEST(LinearInterpolatorTest, ReadsAndMergesPrecursorsWithTolerance)
{
  namespace fs = std::filesystem;

  const auto xs_file = fs::path(OPENSN_TEST_ROOT) / "assets" / "xs" / "simple_fissile.xs";
  const auto xs = MultiGroupXS::LoadFromOpenSn(xs_file.string());
  XSProperties properties(xs);

  ASSERT_EQ(properties.precursors.size(), 1u);
  const auto& [decay_constant, emission_spectrum] = *properties.precursors.begin();
  EXPECT_DOUBLE_EQ(decay_constant, 0.1);
  ASSERT_EQ(emission_spectrum.size(), 1u);
  EXPECT_DOUBLE_EQ(emission_spectrum.front(), 1.0);

  XSProperties nearly_identical;
  const double perturbation = 0.5 * 1e-12;
  nearly_identical.precursors.emplace(0.1 * (1.0 + perturbation),
                                      std::vector<double>{1.0 * (1.0 + perturbation)});
  properties.MergePrecursors(nearly_identical);
  EXPECT_EQ(properties.precursors.size(), 1u);

  XSProperties distinct;
  distinct.precursors.emplace(0.1 * (1.0 + 2.0 * 1e-12), std::vector<double>{1.0});
  properties.MergePrecursors(distinct);
  EXPECT_EQ(properties.precursors.size(), 2u);

  auto with_different_precursors = properties;
  with_different_precursors.precursors.emplace(1.0, std::vector<double>{1.0});
  EXPECT_TRUE(properties == with_different_precursors);
}

TEST(LinearInterpolatorTest, EvaluateInterpolatesMergedPrecursorYields)
{
  const CartesianGrid grid({{0.0, 1.0}});

  auto make_xs = [](double decay_constant, double fractional_yield)
  {
    auto xs = std::make_shared<MultiGroupXS>(1, 0, 1, true);
    xs->GetSigmaTotal() = {1.0};
    xs->GetPrecursors().push_back({decay_constant, fractional_yield, {1.0}});
    return xs;
  };

  const std::vector<std::shared_ptr<MultiGroupXS>> xs = {make_xs(0.1, 0.25), make_xs(0.2, 0.75)};
  LinearInterpolator interpolator(grid, xs, XSType::Total | XSType::Precursor);

  std::vector<double> state_point = {0.5};
  const auto result = interpolator.Evaluate(state_point);

  ASSERT_EQ(result.GetNumPrecursors(), 2u);
  ASSERT_EQ(result.GetPrecursors().size(), 2u);

  const auto& first = result.GetPrecursors()[0];
  EXPECT_DOUBLE_EQ(first.decay_constant, 0.1);
  EXPECT_DOUBLE_EQ(first.emission_spectrum[0], 1.0);
  EXPECT_NEAR(first.fractional_yield, 0.125, 1.0e-12);

  const auto& second = result.GetPrecursors()[1];
  EXPECT_DOUBLE_EQ(second.decay_constant, 0.2);
  EXPECT_DOUBLE_EQ(second.emission_spectrum[0], 1.0);
  EXPECT_NEAR(second.fractional_yield, 0.375, 1.0e-12);
}
