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
