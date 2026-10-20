// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace opensn;

namespace
{

using Matrix = std::vector<std::vector<double>>;

// Two groups, two scattering moments and fission, with values that are exact in binary. The
// transfer matrices (row: to-group, column: from-group) are S0 = [[0.5, 0], [0.25, 0.5]] and
// S1 = [[0, 0], [0.125, 0]], and the production matrix is F = chi nu_sigma_f^T with
// chi = (0.75, 0.25) and nu_sigma_f = (0.25 * 2, 0.5 * 2.5) = (0.5, 1.25).
const std::vector<Matrix> TRANSFER = {{{0.5, 0.0}, {0.25, 0.5}}, {{0.0, 0.0}, {0.125, 0.0}}};
const Matrix PRODUCTION = {{0.375, 0.9375}, {0.125, 0.3125}};

MultiGroupXS
ReadTestXS()
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_xs_adjoint_mode_test_r" + std::to_string(mpi_comm.rank()) + "_" +
                     std::to_string(now) + ".xs");
  {
    std::ofstream out(path);
    out << "NUM_GROUPS 2\n"
           "NUM_MOMENTS 2\n"
           "SIGMA_T_BEGIN\n0 1.0\n1 1.0\nSIGMA_T_END\n"
           "SIGMA_A_BEGIN\n0 0.25\n1 0.5\nSIGMA_A_END\n"
           "SIGMA_F_BEGIN\n0 0.25\n1 0.5\nSIGMA_F_END\n"
           "NU_BEGIN\n0 2.0\n1 2.5\nNU_END\n"
           "CHI_BEGIN\n0 0.75\n1 0.25\nCHI_END\n"
           "TRANSFER_MOMENTS_BEGIN\n"
           "M_GFROM_GTO_VAL 0 0 0 0.5\n"
           "M_GFROM_GTO_VAL 0 0 1 0.25\n"
           "M_GFROM_GTO_VAL 0 1 1 0.5\n"
           "M_GFROM_GTO_VAL 1 0 1 0.125\n"
           "TRANSFER_MOMENTS_END\n";
  }
  auto xs = MultiGroupXS::LoadFromOpenSn(path.string());
  std::filesystem::remove(path);
  return xs;
}

std::shared_ptr<MultiGroupXS>
LoadTestXS()
{
  return std::make_shared<MultiGroupXS>(ReadTestXS());
}

Matrix
Dense(const SparseMatrix& sparse)
{
  Matrix dense(2, std::vector<double>(2, 0.0));
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      dense[i][j] = sparse.GetValueIJ(i, j);
  return dense;
}

Matrix
Scaled(const Matrix& a, const double factor, const bool transpose)
{
  Matrix b(2, std::vector<double>(2, 0.0));
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      b[i][j] = factor * (transpose ? a[j][i] : a[i][j]);
  return b;
}

void
ExpectMatrixEq(const Matrix& actual, const Matrix& expected, const std::string& what)
{
  ASSERT_EQ(actual.size(), expected.size()) << what;
  for (size_t i = 0; i < expected.size(); ++i)
  {
    ASSERT_EQ(actual[i].size(), expected[i].size()) << what;
    for (size_t j = 0; j < expected[i].size(); ++j)
      EXPECT_DOUBLE_EQ(actual[i][j], expected[i][j]) << what << " (" << i << ", " << j << ")";
  }
}

/// Checks the data returned in the set mode against the forward data scaled by factor.
void
ExpectModeData(const MultiGroupXS& xs, const bool adjoint, const double factor)
{
  const std::string mode = adjoint ? "adjoint" : "forward";
  ASSERT_EQ(xs.GetAdjointMode(), adjoint);
  ASSERT_EQ(xs.GetTransferMatrices().size(), TRANSFER.size()) << mode;
  for (unsigned int ell = 0; ell < TRANSFER.size(); ++ell)
  {
    const auto expected = Scaled(TRANSFER[ell], factor, adjoint);
    const auto label = mode + " mode transfer moment " + std::to_string(ell);
    ExpectMatrixEq(Dense(xs.GetTransferMatrix(ell)), expected, label);
    ExpectMatrixEq(Dense(xs.GetTransferMatrices()[ell]), expected, label + " (vector)");
  }
  ExpectMatrixEq(
    xs.GetProductionMatrix(), Scaled(PRODUCTION, factor, adjoint), mode + " mode production");
}

/// Checks the data returned by the getters that take the mode, in both modes, whatever mode is set.
void
ExpectExplicitModeData(const MultiGroupXS& xs, const double factor, const std::string& when)
{
  for (const bool adjoint : {false, true})
  {
    const std::string mode = std::string(adjoint ? "adjoint" : "forward") + " data " + when;
    ASSERT_EQ(xs.GetTransferMatrices(adjoint).size(), TRANSFER.size()) << mode;
    for (unsigned int ell = 0; ell < TRANSFER.size(); ++ell)
    {
      const auto expected = Scaled(TRANSFER[ell], factor, adjoint);
      const auto label = mode + ", transfer moment " + std::to_string(ell);
      ExpectMatrixEq(Dense(xs.GetTransferMatrix(ell, adjoint)), expected, label);
      ExpectMatrixEq(Dense(xs.GetTransferMatrices(adjoint)[ell]), expected, label + " (vector)");
    }
    ExpectMatrixEq(
      xs.GetProductionMatrix(adjoint), Scaled(PRODUCTION, factor, adjoint), mode + ", production");
  }
}

} // namespace

TEST(MultiGroupXSAdjointModeTest, GettersFollowTheSetMode)
{
  auto xs = LoadTestXS();
  ExpectModeData(*xs, false, 1.0);
  xs->SetAdjointMode(true);
  ExpectModeData(*xs, true, 1.0);
  xs->SetAdjointMode(false);
  ExpectModeData(*xs, false, 1.0);
  xs->SetAdjointMode(true);
  ExpectModeData(*xs, true, 1.0);
}

TEST(MultiGroupXSAdjointModeTest, EnteringAdjointModeTransposesDataScaledInForwardMode)
{
  auto xs = LoadTestXS();
  xs->SetAdjointMode(true);
  xs->SetAdjointMode(false);
  xs->Scale(2.0);
  ExpectModeData(*xs, false, 2.0);
  xs->SetAdjointMode(true);
  ExpectModeData(*xs, true, 2.0);

  // Scaling in adjoint mode keeps the transposes current.
  xs->Scale(3.0);
  ExpectModeData(*xs, true, 3.0);
  xs->SetAdjointMode(false);
  ExpectModeData(*xs, false, 3.0);
}

TEST(MultiGroupXSAdjointModeTest, GettersWithModeFollowScaleInEitherMode)
{
  // Once adjoint mode has been set, the transposes must follow every Scale, also while the
  // getters without a mode argument are in forward mode.
  auto xs = LoadTestXS();
  xs->SetAdjointMode(true);
  xs->SetAdjointMode(false);
  for (const double factor : {2.0, 3.0})
  {
    xs->Scale(factor);
    ExpectModeData(*xs, false, factor);
    const auto when =
      "after Scale(" + std::to_string(static_cast<int>(factor)) + ") in forward mode";
    ExpectExplicitModeData(*xs, factor, when);
  }

  xs->SetAdjointMode(true);
  ExpectExplicitModeData(*xs, 3.0, "after entering adjoint mode");
  xs->Scale(4.0);
  ExpectExplicitModeData(*xs, 4.0, "after Scale(4) in adjoint mode");
  xs->SetAdjointMode(false);
  ExpectExplicitModeData(*xs, 4.0, "after leaving adjoint mode");
}

TEST(MultiGroupXSAdjointModeTest, ReplaceContentsKeepsAdjointModeAndTransposesTheNewData)
{
  // A problem in adjoint mode that holds this object reads the transposes. Replacing the contents
  // with new data, as the Python methods that reload cross sections in place do, must keep the
  // adjoint mode and give the transposes of the new data, also after the new data are scaled.
  auto xs = LoadTestXS();
  xs->Scale(3.0);
  xs->SetAdjointMode(true);
  ExpectModeData(*xs, true, 3.0);

  xs->ReplaceContents(ReadTestXS());
  ExpectModeData(*xs, true, 1.0);
  ExpectExplicitModeData(*xs, 1.0, "after ReplaceContents");

  xs->Scale(2.0);
  ExpectModeData(*xs, true, 2.0);
  ExpectExplicitModeData(*xs, 2.0, "after Scale(2) of the new data");
}

TEST(MultiGroupXSAdjointModeTest, ReplaceContentsWithACopyKeepsTheModesOfBothObjects)
{
  // An object in adjoint mode takes the data of a copy of an object in forward mode.
  auto xs = LoadTestXS();
  xs->Scale(3.0);
  xs->SetAdjointMode(true);
  auto forward = LoadTestXS();
  forward->Scale(0.5);
  xs->ReplaceContents(MultiGroupXS(*forward));
  ExpectModeData(*xs, true, 0.5);
  ExpectExplicitModeData(*xs, 0.5, "after ReplaceContents with a forward-mode copy");
  ExpectModeData(*forward, false, 0.5);

  // An object in forward mode takes the data of a copy of an object in adjoint mode.
  auto target = LoadTestXS();
  auto adjoint_source = LoadTestXS();
  adjoint_source->Scale(0.5);
  adjoint_source->SetAdjointMode(true);
  target->ReplaceContents(MultiGroupXS(*adjoint_source));
  ExpectModeData(*target, false, 0.5);
  ExpectModeData(*adjoint_source, true, 0.5);

  // An object in adjoint mode takes the data of a copy of another object in adjoint mode. The copy
  // brings its own transposes, which must be replaced, not added to.
  auto adjoint_target = LoadTestXS();
  adjoint_target->Scale(3.0);
  adjoint_target->SetAdjointMode(true);
  adjoint_target->ReplaceContents(MultiGroupXS(*adjoint_source));
  ExpectModeData(*adjoint_target, true, 0.5);
  ExpectExplicitModeData(*adjoint_target, 0.5, "after ReplaceContents with an adjoint-mode copy");
  ExpectModeData(*adjoint_source, true, 0.5);
}

TEST(MultiGroupXSAdjointModeTest, ReplaceContentsInForwardModeRebuildsExistingTransposes)
{
  // Transposes built in adjoint mode stay current after the object returns to forward mode and
  // its contents are replaced.
  auto xs = LoadTestXS();
  xs->Scale(3.0);
  xs->SetAdjointMode(true);
  xs->SetAdjointMode(false);

  xs->ReplaceContents(ReadTestXS());
  ExpectModeData(*xs, false, 1.0);
  ExpectExplicitModeData(*xs, 1.0, "after ReplaceContents in forward mode");

  xs->Scale(2.0);
  ExpectExplicitModeData(*xs, 2.0, "after Scale(2) of the new data in forward mode");
}

TEST(MultiGroupXSAdjointModeTest, CombineUsesForwardData)
{
  auto xs_adjoint = LoadTestXS();
  auto xs_forward = LoadTestXS();
  xs_adjoint->SetAdjointMode(true);

  for (const auto& combined : {MultiGroupXS::Combine({{xs_adjoint, 2.0}}),
                               MultiGroupXS::Combine({{xs_adjoint, 1.0}, {xs_forward, 1.0}}),
                               MultiGroupXS::Combine({{xs_forward, 1.0}, {xs_adjoint, 1.0}})})
    ExpectModeData(combined, false, 2.0);
}

TEST(MultiGroupXSAdjointModeTest, CollisionEnergyLossUsesForwardTransfer)
{
  // coefficient(g) = E_g sigma_t,g - sum_g' E_g' sigma_s,0(g -> g') with group centers (2, 0.5):
  // (2 - 2 * 0.5 - 0.5 * 0.25, 0.5 - 0.5 * 0.5) = (0.875, 0.25).
  auto xs = LoadTestXS();
  const EnergyGroupStructure energy{{2.0, 0.5}, {1.0, 1.0}};
  for (const bool adjoint : {false, true})
  {
    xs->SetAdjointMode(adjoint);
    const auto coefficients = xs->ComputeCollisionEnergyLossCoefficients(energy);
    ASSERT_EQ(coefficients.size(), 2U);
    EXPECT_DOUBLE_EQ(coefficients[0], 0.875) << "adjoint mode " << adjoint;
    EXPECT_DOUBLE_EQ(coefficients[1], 0.25) << "adjoint mode " << adjoint;
  }
}
