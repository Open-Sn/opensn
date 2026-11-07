// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/acceleration.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

using namespace opensn;

namespace
{

/**
 * Loads a cross section with the given total cross sections and isotropic transfer matrix
 * transfer[g_to][g_from] from a temporary OpenSn cross-section file.
 */
MultiGroupXS
MakeXS(const std::vector<double>& sigma_t, const std::vector<std::vector<double>>& transfer)
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_two_grid_spectrum_r" + std::to_string(mpi_comm.rank()) + "_" +
                     std::to_string(now) + ".xs");
  {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "NUM_GROUPS " << sigma_t.size() << "\nNUM_MOMENTS 1\n";
    out << "SIGMA_T_BEGIN\n";
    for (std::size_t g = 0; g < sigma_t.size(); ++g)
      out << g << " " << sigma_t[g] << "\n";
    out << "SIGMA_T_END\n";
    out << "TRANSFER_MOMENTS_BEGIN\n";
    for (std::size_t g_to = 0; g_to < transfer.size(); ++g_to)
      for (std::size_t g_from = 0; g_from < transfer[g_to].size(); ++g_from)
        if (transfer[g_to][g_from] != 0.0)
          out << "M_GFROM_GTO_VAL 0 " << g_from << " " << g_to << " " << transfer[g_to][g_from]
              << "\n";
    out << "TRANSFER_MOMENTS_END\n";
  }
  auto xs = MultiGroupXS::LoadFromOpenSn(path.string());
  std::filesystem::remove(path);
  return xs;
}

void
ExpectSpectrum(const std::vector<double>& spectrum, const std::vector<double>& expected)
{
  ASSERT_EQ(spectrum.size(), expected.size());
  for (std::size_t g = 0; g < expected.size(); ++g)
    EXPECT_NEAR(spectrum[g], expected[g], 1.0e-10) << "group " << g;
}

} // namespace

// Upscatter only, with converged within-group scattering (JFULL): the iteration matrix is strictly
// upper triangular, hence nilpotent, so there is no error mode to remove and the spectrum is flat.
TEST(TwoGridSpectrumTest, UpscatterOnlyIsFlat)
{
  const auto xs = MakeXS({1.0, 1.0, 1.0}, {{0.5, 0.3, 0.2}, {0.0, 0.5, 0.3}, {0.0, 0.0, 0.5}});
  const auto info = MakeTwoGridCollapsedInfo(xs, EnergyCollapseScheme::JFULL, 0, 2);
  ExpectSpectrum(info.spectrum, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
}

// Downscatter only (JFULL): strictly lower triangular, nilpotent, flat spectrum.
TEST(TwoGridSpectrumTest, DownscatterOnlyIsFlat)
{
  const auto xs = MakeXS({1.0, 1.0, 1.0}, {{0.5, 0.0, 0.0}, {0.3, 0.5, 0.0}, {0.2, 0.3, 0.5}});
  const auto info = MakeTwoGridCollapsedInfo(xs, EnergyCollapseScheme::JFULL, 0, 2);
  ExpectSpectrum(info.spectrum, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
}

// Two groups coupled by down- and upscatter (JFULL). With sigma_t = (1, 2) and
// S = [[0.5, 0.09], [0.45, 1.9]], the iteration matrix (diag(sigma_t) - diag(S))^-1 offdiag(S) is
// [[0, 0.18], [4.5, 0]] with eigenvalues +/-0.9. The Perron eigenvector satisfies
// 0.18 xi_1 = 0.9 xi_0, so xi = (1/6, 5/6).
TEST(TwoGridSpectrumTest, TwoGroupPerronVector)
{
  const auto xs = MakeXS({1.0, 2.0}, {{0.5, 0.09}, {0.45, 1.9}});
  const auto info = MakeTwoGridCollapsedInfo(xs, EnergyCollapseScheme::JFULL, 0, 1);
  ExpectSpectrum(info.spectrum, {1.0 / 6.0, 5.0 / 6.0});
}

namespace
{

/**
 * Checks the defining identity of the local two-grid correction, L (Sigma_t - S) = h(C) with
 * h(C) = (rho I - C) / (2 rho), for the Jacobi iteration matrix
 * C = (diag(sigma_t) - diag(S))^-1 offdiag(S) and its Perron eigenvalue rho, and that the
 * diffusion-source weights c = (I - (Sigma_t - S) L)^T 1 leave the balance of the Perron mode xi
 * unchanged: c^T (Sigma_t - S) xi = 1^T (Sigma_t - S) xi.
 */
void
ExpectLocalCorrection(const TwoGridCollapsedInfo& info,
                      const std::vector<double>& sigma_t,
                      const std::vector<std::vector<double>>& transfer,
                      const double rho)
{
  const auto n = sigma_t.size();
  ASSERT_EQ(info.local_correction.size(), n * n);
  ASSERT_EQ(info.source_weights.size(), n);
  for (std::size_t g = 0; g < n; ++g)
    for (std::size_t gp = 0; gp < n; ++gp)
    {
      double L_AmB = 0.0;
      for (std::size_t k = 0; k < n; ++k)
        L_AmB +=
          info.local_correction[g * n + k] * ((k == gp ? sigma_t[k] : 0.0) - transfer[k][gp]);
      const double C = g == gp ? 0.0 : transfer[g][gp] / (sigma_t[g] - transfer[g][g]);
      const double h = ((g == gp ? rho : 0.0) - C) / (2.0 * rho);
      EXPECT_NEAR(L_AmB, h, 1.0e-10) << "entry (" << g << ", " << gp << ")";
    }

  double c_balance = 0.0;
  double balance = 0.0;
  for (std::size_t g = 0; g < n; ++g)
  {
    double AmB_xi = sigma_t[g] * info.spectrum[g];
    for (std::size_t gp = 0; gp < n; ++gp)
      AmB_xi -= transfer[g][gp] * info.spectrum[gp];
    c_balance += info.source_weights[g] * AmB_xi;
    balance += AmB_xi;
  }
  EXPECT_NEAR(c_balance, balance, 1.0e-12);
}

} // namespace

// Two groups with Jacobi eigenvalues +/-0.9 (TwoGroupPerronVector). The local correction removes
// the -0.9 mode, xi_- = (1, -5), exactly: L (Sigma_t - S) xi_- = xi_-, and leaves the Perron mode
// xi_+ = (1, 5) to the diffusion correction: L (Sigma_t - S) xi_+ = 0.
TEST(TwoGridSpectrumTest, TwoGroupAlternatingModeIsRemovedLocally)
{
  const std::vector<double> sigma_t = {1.0, 2.0};
  const std::vector<std::vector<double>> transfer = {{0.5, 0.09}, {0.45, 1.9}};
  const auto info =
    MakeTwoGridCollapsedInfo(MakeXS(sigma_t, transfer), EnergyCollapseScheme::JFULL, 0, 1);
  ExpectLocalCorrection(info, sigma_t, transfer, 0.9);

  // (Sigma_t - S) xi_- = (0.95, -0.95) and (Sigma_t - S) xi_+ = (0.05, 0.05). h(C) cancels the
  // Perron mode, amplified by 1 / (1 - rho) = 10, only to the accuracy of rho, about 1e-13.
  const auto& L = info.local_correction;
  EXPECT_NEAR(L[0] * 0.95 - L[1] * 0.95, 1.0, 1.0e-10);
  EXPECT_NEAR(L[2] * 0.95 - L[3] * 0.95, -5.0, 1.0e-10);
  EXPECT_NEAR(L[0] * 0.05 + L[1] * 0.05, 0.0, 1.0e-10);
  EXPECT_NEAR(L[2] * 0.05 + L[3] * 0.05, 0.0, 1.0e-10);
}

// Three groups with equal coupling: C = 0.4 (J - I) has the Perron eigenvalue 0.8 (flat spectrum)
// and the double eigenvalue -0.4, whose whole eigenspace receives the weight h(-0.4) = 0.75.
TEST(TwoGridSpectrumTest, ThreeGroupDoubleNegativeEigenvalue)
{
  const std::vector<double> sigma_t = {1.0, 1.0, 1.0};
  const std::vector<std::vector<double>> transfer = {
    {0.5, 0.2, 0.2}, {0.2, 0.5, 0.2}, {0.2, 0.2, 0.5}};
  const auto info =
    MakeTwoGridCollapsedInfo(MakeXS(sigma_t, transfer), EnergyCollapseScheme::JFULL, 0, 2);
  ExpectSpectrum(info.spectrum, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
  ExpectLocalCorrection(info, sigma_t, transfer, 0.8);
}

// Two groups without absorption: rho = 1 and Sigma_t - S is singular, so there is no local
// correction and the two-grid correction is that of the spectrum alone.
TEST(TwoGridSpectrumTest, NoAbsorptionHasNoLocalCorrection)
{
  const auto xs = MakeXS({1.0, 2.0}, {{0.5, 0.1}, {0.5, 1.9}});
  const auto info = MakeTwoGridCollapsedInfo(xs, EnergyCollapseScheme::JFULL, 0, 1);
  ExpectSpectrum(info.spectrum, {1.0 / 6.0, 5.0 / 6.0});
  EXPECT_TRUE(info.local_correction.empty());
  ASSERT_EQ(info.source_weights.size(), 2U);
  for (const auto w : info.source_weights)
    EXPECT_EQ(w, 1.0);
}

// Upscatter only: the iteration matrix is nilpotent, there is no Perron mode, and neither the
// local correction nor the source weights change the two-grid correction.
TEST(TwoGridSpectrumTest, NilpotentHasNoLocalCorrection)
{
  const auto xs = MakeXS({1.0, 1.0, 1.0}, {{0.5, 0.3, 0.2}, {0.0, 0.5, 0.3}, {0.0, 0.0, 0.5}});
  const auto info = MakeTwoGridCollapsedInfo(xs, EnergyCollapseScheme::JFULL, 0, 2);
  EXPECT_TRUE(info.local_correction.empty());
  ASSERT_EQ(info.source_weights.size(), 3U);
  for (const auto w : info.source_weights)
    EXPECT_EQ(w, 1.0);
}
