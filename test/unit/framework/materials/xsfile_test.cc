// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/xsfile.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <array>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <memory>
#include <vector>
#include <stdexcept>

using namespace opensn;

namespace
{

std::filesystem::path
MakeTempPath(const std::string& stem, const std::string& extension)
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  return std::filesystem::temp_directory_path() /
         (stem + "_r" + std::to_string(mpi_comm.rank()) + "_" + std::to_string(now) + extension);
}

void
WriteU32(std::ostream& out, const std::uint32_t value)
{
  std::array<char, sizeof(value)> bytes{};
  std::memcpy(bytes.data(), &value, sizeof(value));
  out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

/// Writes a two-group OpenSn cross-section file. Empty bounds omit GROUP_STRUCTURE. The extra
/// text is appended verbatim.
std::string
MakeTempXSFile(const std::vector<double>& bounds = {1.0, 0.1, 0.01}, const std::string& extra = "")
{
  const auto path = MakeTempPath("opensn_xsfile_test", ".xs");

  std::ofstream out(path);
  out << std::setprecision(17);
  out << "NUM_GROUPS 2\n";
  if (not bounds.empty())
  {
    out << "GROUP_STRUCTURE_BEGIN\n";
    for (const double bound : bounds)
      out << bound << "\n";
    out << "GROUP_STRUCTURE_END\n";
  }
  out << "SIGMA_T_BEGIN\n";
  out << "0 1.0\n";
  out << "1 2.0\n";
  out << "SIGMA_T_END\n";
  out << "SIGMA_A_BEGIN\n";
  out << "0 0.5\n";
  out << "1 1.5\n";
  out << "SIGMA_A_END\n";
  out << extra;
  out.close();

  return path.string();
}

} // namespace

TEST(XSFileTest, ReadMinimalFile)
{
  const std::string fname = MakeTempXSFile();
  XSFile xs(fname);
  xs.Read();

  EXPECT_EQ(xs.num_groups_, 2U);
  ASSERT_EQ(xs.e_bounds_.size(), 3U);
  EXPECT_NEAR(xs.e_bounds_[0], 1.0, 1e-12);
  EXPECT_NEAR(xs.e_bounds_[1], 0.1, 1e-12);
  EXPECT_NEAR(xs.e_bounds_[2], 0.01, 1e-12);

  ASSERT_EQ(xs.sigma_t_.size(), 2U);
  EXPECT_NEAR(xs.sigma_t_[0], 1.0, 1e-12);
  EXPECT_NEAR(xs.sigma_t_[1], 2.0, 1e-12);

  ASSERT_EQ(xs.sigma_a_.size(), 2U);
  EXPECT_NEAR(xs.sigma_a_[0], 0.5, 1e-12);
  EXPECT_NEAR(xs.sigma_a_[1], 1.5, 1e-12);

  std::filesystem::remove(fname);
}

TEST(MultiGroupXSTest, CombineAllowsInputsWithoutEnergyBounds)
{
  const auto bounded_file = MakeTempXSFile();
  const auto unbounded_file = MakeTempXSFile({});
  const auto other_bounds_file = MakeTempXSFile({2.0, 0.1, 0.01});
  auto bounded = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(bounded_file));
  auto unbounded = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(unbounded_file));
  auto other_bounds =
    std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(other_bounds_file));
  std::filesystem::remove(bounded_file);
  std::filesystem::remove(unbounded_file);
  std::filesystem::remove(other_bounds_file);
  ASSERT_TRUE(bounded->HasEnergyGroupBounds());
  ASSERT_FALSE(unbounded->HasEnergyGroupBounds());

  // Bounds come from whichever input supplies them, regardless of order.
  for (const auto& combined : {MultiGroupXS::Combine({{bounded, 1.0}, {unbounded, 1.0}}),
                               MultiGroupXS::Combine({{unbounded, 1.0}, {bounded, 1.0}})})
  {
    ASSERT_TRUE(combined.HasEnergyGroupBounds());
    for (unsigned int g = 0; g < 2; ++g)
      EXPECT_EQ(combined.GetEnergyGroupBounds(g), bounded->GetEnergyGroupBounds(g));
    EXPECT_DOUBLE_EQ(combined.GetSigmaTotal()[1], 4.0);
  }

  const auto unbounded_only = MultiGroupXS::Combine({{unbounded, 1.0}, {unbounded, 1.0}});
  EXPECT_FALSE(unbounded_only.HasEnergyGroupBounds());

  EXPECT_THROW(MultiGroupXS::Combine({{bounded, 1.0}, {unbounded, 1.0}, {other_bounds, 1.0}}),
               std::logic_error);
}

TEST(MultiGroupXSTest, EnergyBoundsMatchWithinRoundoff)
{
  // Bounds that differ only by round-off must be accepted; a real mismatch must not.
  const auto bounded_file = MakeTempXSFile();
  const auto nearby_file = MakeTempXSFile({1.0 * (1.0 + 1.0e-9), 0.1 * (1.0 - 1.0e-9), 0.01});
  const auto different_file = MakeTempXSFile({1.0, 0.1 * (1.0 + 1.0e-3), 0.01});
  auto bounded = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(bounded_file));
  auto nearby = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(nearby_file));
  auto different = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(different_file));
  std::filesystem::remove(bounded_file);
  std::filesystem::remove(nearby_file);
  std::filesystem::remove(different_file);
  ASSERT_NE(bounded->GetEnergyGroupBounds(0), nearby->GetEnergyGroupBounds(0));

  EXPECT_NO_THROW(MultiGroupXS::Combine({{bounded, 1.0}, {nearby, 1.0}}));
  EXPECT_THROW(MultiGroupXS::Combine({{bounded, 1.0}, {different, 1.0}}), std::logic_error);

  EXPECT_NO_THROW(MultiGroupXS::ResolveEnergyGroupStructure({{0, bounded}, {1, nearby}}, 2));
  EXPECT_THROW(MultiGroupXS::ResolveEnergyGroupStructure({{0, bounded}, {1, different}}, 2),
               std::invalid_argument);
}

TEST(MultiGroupXSTest, AdjointModeWithoutTransferMatrices)
{
  // A file without a TRANSFER_MOMENTS block describes a material with no scattering, so it has
  // no transfer matrices. Its adjoint has no scattering either.
  const auto fname = MakeTempXSFile();
  auto xs = MultiGroupXS::LoadFromOpenSn(fname);
  std::filesystem::remove(fname);
  ASSERT_TRUE(xs.GetTransferMatrices().empty());

  for (const bool adjoint : {true, false, true})
  {
    xs.SetAdjointMode(adjoint);
    EXPECT_EQ(xs.GetAdjointMode(), adjoint);
    for (const auto& S_ell : xs.GetTransferMatrices())
      for (unsigned int g = 0; g < xs.GetNumGroups(); ++g)
        for (const auto& entry : S_ell.Row(g))
          EXPECT_EQ(entry.value, 0.0);
    ASSERT_EQ(xs.GetSigmaTotal().size(), 2U);
    EXPECT_EQ(xs.GetSigmaTotal()[0], 1.0);
    EXPECT_EQ(xs.GetSigmaTotal()[1], 2.0);
  }
}

TEST(MultiGroupXSTest, AdjointProductionWithoutTransferMatrices)
{
  // A fissionable material without scattering has no transfer matrices, but adjoint mode must
  // still transpose its production matrix, and scaling must scale the transpose. NUM_MOMENTS 2
  // without a TRANSFER_MOMENTS block gives scattering order 1 and still no transfer matrices.
  const auto fname = MakeTempXSFile({1.0, 0.1, 0.01},
                                    "NUM_MOMENTS 2\n"
                                    "SIGMA_F_BEGIN\n0 0.25\n1 0.5\nSIGMA_F_END\n"
                                    "NU_BEGIN\n0 2.5\n1 2.0\nNU_END\n"
                                    "CHI_BEGIN\n0 0.75\n1 0.25\nCHI_END\n");
  auto xs = MultiGroupXS::LoadFromOpenSn(fname);
  std::filesystem::remove(fname);
  ASSERT_TRUE(xs.IsFissionable());
  ASSERT_TRUE(xs.GetTransferMatrices().empty());

  // F[g][g'] = chi[g] * nu_sigma_f[g'], with chi = (0.75, 0.25) and nu_sigma_f = (0.625, 1.0).
  // These values and their doubles are exact in binary.
  using Matrix = std::vector<std::vector<double>>;
  const Matrix F = {{0.46875, 0.75}, {0.15625, 0.25}};
  const Matrix F_transpose = {{0.46875, 0.15625}, {0.75, 0.25}};
  const Matrix twice_F = {{0.9375, 1.5}, {0.3125, 0.5}};
  const Matrix twice_F_transpose = {{0.9375, 0.3125}, {1.5, 0.5}};
  EXPECT_EQ(xs.GetProductionMatrix(), F);

  xs.SetAdjointMode(true);
  EXPECT_TRUE(xs.GetTransferMatrices().empty());
  EXPECT_EQ(xs.GetProductionMatrix(), F_transpose);

  xs.Scale(2.0);
  EXPECT_TRUE(xs.GetAdjointMode());
  EXPECT_TRUE(xs.GetTransferMatrices().empty());
  EXPECT_EQ(xs.GetNuSigmaF(), (std::vector<double>{1.25, 2.0}));
  EXPECT_EQ(xs.GetProductionMatrix(), twice_F_transpose);

  xs.SetAdjointMode(false);
  EXPECT_EQ(xs.GetProductionMatrix(), twice_F);
  xs.SetAdjointMode(true);
  EXPECT_EQ(xs.GetProductionMatrix(), twice_F_transpose);
}

TEST(MultiGroupXSTest, AdjointModeTransposesEveryTransferMoment)
{
  // Adjoint mode transposes every transfer matrix that the file defines, here P0 and P1, and
  // scaling in adjoint mode scales them.
  const auto fname = MakeTempXSFile({1.0, 0.1, 0.01},
                                    "NUM_MOMENTS 2\n"
                                    "TRANSFER_MOMENTS_BEGIN\n"
                                    "M_GFROM_GTO_VAL 0 0 0 0.25\n"
                                    "M_GFROM_GTO_VAL 0 0 1 0.125\n"
                                    "M_GFROM_GTO_VAL 0 1 1 0.375\n"
                                    "M_GFROM_GTO_VAL 1 0 0 0.0625\n"
                                    "M_GFROM_GTO_VAL 1 0 1 0.03125\n"
                                    "M_GFROM_GTO_VAL 1 1 1 0.09375\n"
                                    "TRANSFER_MOMENTS_END\n");
  auto xs = MultiGroupXS::LoadFromOpenSn(fname);
  std::filesystem::remove(fname);

  // S[ell][g][g'] is the moment-ell transfer from group g' to group g.
  const std::vector<std::vector<std::vector<double>>> S = {{{0.25, 0.0}, {0.125, 0.375}},
                                                           {{0.0625, 0.0}, {0.03125, 0.09375}}};
  const auto check = [&](const bool adjoint, const double factor)
  {
    ASSERT_EQ(xs.GetTransferMatrices().size(), 2U);
    for (unsigned int ell = 0; ell < 2; ++ell)
      for (unsigned int g = 0; g < 2; ++g)
        for (unsigned int gp = 0; gp < 2; ++gp)
          EXPECT_EQ(xs.GetTransferMatrix(ell).GetValueIJ(g, gp),
                    factor * (adjoint ? S[ell][gp][g] : S[ell][g][gp]));
  };

  check(false, 1.0);
  xs.SetAdjointMode(true);
  check(true, 1.0);
  xs.Scale(2.0);
  check(true, 2.0);
  xs.SetAdjointMode(false);
  check(false, 2.0);
}

TEST(MultiGroupXSTest, AdjointModeWithNumMomentsAfterTransferMoments)
{
  // The reader sizes the TRANSFER_MOMENTS block by the NUM_MOMENTS read before it, or one moment
  // if none was. With NUM_MOMENTS 2 after the block, the file yields one transfer matrix and
  // scattering order 1. Adjoint mode must transpose the transfer matrix that exists, and scaling
  // in adjoint mode must scale it.
  const auto fname = MakeTempXSFile({1.0, 0.1, 0.01},
                                    "TRANSFER_MOMENTS_BEGIN\n"
                                    "M_GFROM_GTO_VAL 0 0 0 0.25\n"
                                    "M_GFROM_GTO_VAL 0 0 1 0.125\n"
                                    "M_GFROM_GTO_VAL 0 1 1 0.375\n"
                                    "TRANSFER_MOMENTS_END\n"
                                    "NUM_MOMENTS 2\n");
  auto xs = MultiGroupXS::LoadFromOpenSn(fname);
  std::filesystem::remove(fname);
  ASSERT_EQ(xs.GetScatteringOrder(), 1U);

  // S[g][g'] is the P0 transfer from group g' to group g.
  const std::vector<std::vector<double>> S = {{0.25, 0.0}, {0.125, 0.375}};
  const auto check = [&](const bool adjoint, const double factor)
  {
    ASSERT_EQ(xs.GetTransferMatrices().size(), 1U);
    for (unsigned int g = 0; g < 2; ++g)
      for (unsigned int gp = 0; gp < 2; ++gp)
        EXPECT_EQ(xs.GetTransferMatrix(0).GetValueIJ(g, gp),
                  factor * (adjoint ? S[gp][g] : S[g][gp]));
  };

  check(false, 1.0);
  xs.SetAdjointMode(true);
  check(true, 1.0);
  xs.Scale(2.0);
  check(true, 2.0);
  xs.SetAdjointMode(false);
  check(false, 2.0);
}

TEST(CEPXS, CoupledEnergyBounds)
{
  // In this 40-electron/40-photon library both species share the upper energy.
  // The first photon group must not start at the electron cutoff (0.01 MeV).
  const auto path = std::string(OPENSN_TEST_ROOT) + "/assets/xs/Cu_40ge_40gp_p15_CEPXS_CSDA.bxslib";
  EXPECT_THROW(MultiGroupXS::LoadFromCEPXS(path, 0, true), std::invalid_argument);
  const auto xs = MultiGroupXS::LoadFromCEPXS(path, 0, true, {"electron", "photon"});
  EXPECT_TRUE(xs.HasEnergyGroupBounds());
  const auto [upper_e, lower_e] = xs.GetEnergyGroupBounds(0);
  const auto [upper_p, lower_p] = xs.GetEnergyGroupBounds(40);
  EXPECT_DOUBLE_EQ(upper_p, upper_e);
  EXPECT_NEAR(lower_p, 0.987468344135375, 1.e-14);
  EXPECT_NEAR(0.5 * (upper_p + lower_p), 1.0, 1.e-8);
  EXPECT_NEAR(xs.GetDeltaE()[40], upper_e - lower_p, 1.e-14);
  EXPECT_GT(upper_e, lower_e);
  EXPECT_THROW(xs.GetEnergyGroupBounds(xs.GetNumGroups()), std::out_of_range);

  const auto& particle_types = xs.GetParticleTypes();
  ASSERT_EQ(particle_types.size(), 80);
  EXPECT_TRUE(std::all_of(particle_types.begin(),
                          particle_types.begin() + 40,
                          [](const auto type) { return type == ParticleType::ELECTRON; }));
  EXPECT_TRUE(std::all_of(particle_types.begin() + 40,
                          particle_types.end(),
                          [](const auto type) { return type == ParticleType::PHOTON; }));

  EXPECT_NO_THROW(MultiGroupXS::LoadFromCEPXS(path, 0, true, {"electron", "photon"}));
  EXPECT_THROW(MultiGroupXS::LoadFromCEPXS(path, 0, true, {"electron"}), std::invalid_argument);
  EXPECT_THROW(MultiGroupXS::LoadFromCEPXS(path, 0, true, {"photon", "electron"}),
               std::invalid_argument);
  EXPECT_THROW(MultiGroupXS::LoadFromCEPXS(path, 0, true, {"electron", "electron"}),
               std::invalid_argument);
  EXPECT_THROW(MultiGroupXS::LoadFromCEPXS(path, 0, true, {"electron", "muon"}),
               std::invalid_argument);
}

TEST(CEPXS, RejectsMalformedFortranRecords)
{
  const auto fixture =
    std::filesystem::path(OPENSN_TEST_ROOT) / "assets/xs/cepxs_synthetic_csda_3g.bxslib";

  const auto partial_header = MakeTempPath("opensn_cepxs_partial_header", ".bxslib");
  std::filesystem::copy_file(fixture, partial_header);
  {
    std::ofstream out(partial_header, std::ios::binary | std::ios::app);
    out.put('\0');
  }
  EXPECT_THROW(MultiGroupXS::LoadFromCEPXS(partial_header.string(), 0, true), std::runtime_error);
  std::filesystem::remove(partial_header);

  const auto oversized_record = MakeTempPath("opensn_cepxs_oversized_record", ".bxslib");
  {
    std::ofstream out(oversized_record, std::ios::binary);
    const std::uint32_t title_size = 1;
    WriteU32(out, title_size);
    out.put('X');
    WriteU32(out, title_size);
    const std::uint32_t oversized_marker = std::numeric_limits<std::uint32_t>::max();
    WriteU32(out, oversized_marker);
  }
  EXPECT_THROW(MultiGroupXS::LoadFromCEPXS(oversized_record.string(), 0, true), std::runtime_error);
  std::filesystem::remove(oversized_record);
}

TEST(CEPXS, ScaleAndCombine)
{
  for (const bool csda : {false, true})
  {
    SCOPED_TRACE(csda);
    const auto original = MultiGroupXS::LoadFromCEPXS(
      std::string(OPENSN_TEST_ROOT) + "/assets/xs/Cu_40ge_40gp_p15_CEPXS_CSDA.bxslib",
      0,
      csda,
      csda ? std::vector<std::string>{"electron", "photon"} : std::vector<std::string>{});
    // Only the CSDA row layout carries stopping power; legacy imports leave it absent.
    EXPECT_EQ(original.GetStoppingPower().empty(), not csda);
    EXPECT_EQ(original.GetByName("stopping_power") == nullptr, not csda);
    ASSERT_TRUE(original.HasCustomXS("cepxs_charge_deposition"));
    ASSERT_TRUE(original.HasCustomXS("cepxs_secondary_production"));
    EXPECT_EQ(original.GetCustomXS("charge_deposition"),
              original.GetCustomXS("cepxs_charge_deposition"));

    const auto CheckScaled = [&](const MultiGroupXS& xs, const double factor)
    {
      const auto CheckVector = [factor](const auto& actual, const auto& base)
      {
        ASSERT_EQ(actual.size(), base.size());
        for (size_t i = 0; i < base.size(); ++i)
          EXPECT_NEAR(
            actual[i], factor * base[i], 1.e-12 * std::max(1.0, std::abs(factor * base[i])));
      };
      CheckVector(xs.GetSigmaTotal(), original.GetSigmaTotal());
      CheckVector(xs.GetSigmaAbsorption(), original.GetSigmaAbsorption());
      CheckVector(xs.GetEnergyDeposition(), original.GetEnergyDeposition());
      CheckVector(xs.GetStoppingPower(), original.GetStoppingPower());
      for (const auto& name : original.GetCustomXSNames())
        CheckVector(xs.GetCustomXS(name), original.GetCustomXS(name));
      ASSERT_EQ(xs.GetTransferMatrices().size(), original.GetTransferMatrices().size());
      for (size_t m = 0; m < original.GetTransferMatrices().size(); ++m)
        for (unsigned int g = 0; g < original.GetNumGroups(); ++g)
        {
          const auto& actual = xs.GetTransferMatrices()[m];
          const auto& base = original.GetTransferMatrices()[m];
          EXPECT_EQ(actual.rowI_indices[g], base.rowI_indices[g]);
          CheckVector(actual.rowI_values[g], base.rowI_values[g]);
        }
      for (unsigned int g = 0; g < original.GetNumGroups(); ++g)
        EXPECT_EQ(xs.GetEnergyGroupBounds(g), original.GetEnergyGroupBounds(g));
      EXPECT_EQ(xs.GetDeltaE(), original.GetDeltaE());
    };

    auto first = std::make_shared<MultiGroupXS>(original);
    first->Scale(2.0);
    CheckScaled(*first, 2.0);
    first->Scale(3.0);
    CheckScaled(*first, 3.0); // Factors are absolute, not cumulative.
    auto second = std::make_shared<MultiGroupXS>(original);
    second->Scale(0.5);
    auto combined = MultiGroupXS::Combine({{first, 0.25}, {second, 0.5}});
    CheckScaled(combined, 1.0); // 0.25 * 3 + 0.5 * 0.5 = 1.
    combined.Scale(2.0);
    CheckScaled(combined, 2.0);
    combined.Scale(1.0);
    CheckScaled(combined, 1.0);
    CheckScaled(*first, 3.0);
    CheckScaled(*second, 0.5);
    first->Scale(0.0);
    CheckScaled(*first, 0.0);
    first->Scale(1.0);
    CheckScaled(*first, 1.0);
  }
}
