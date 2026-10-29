// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/utils/hdf_utils.h"
#include "framework/runtime.h"
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#include <gtest/gtest.h>

using namespace opensn;

namespace
{

constexpr double kTolerance = 1.0e-12;

/// A dataset of an OpenMC cross-section set: row-major values and dimensions.
struct Dataset
{
  std::vector<double> values;
  std::vector<hsize_t> dims;
};

using Datasets = std::map<std::string, std::optional<Dataset>>;

Dataset
Dataset1D(const std::vector<double>& values)
{
  return {values, {values.size()}};
}

Dataset
Dataset2D(const std::vector<double>& values, const hsize_t rows, const hsize_t columns)
{
  return {values, {rows, columns}};
}

bool
WriteStringAttribute(hid_t id, const std::string& name, const std::string& value)
{
  const hid_t type = H5Tcopy(H5T_C_S1);
  if (type < 0)
    return false;
  H5Tset_size(type, value.size());
  H5Tset_strpad(type, H5T_STR_NULLTERM);

  const hid_t dataspace = H5Screate(H5S_SCALAR);
  if (dataspace < 0)
  {
    H5Tclose(type);
    return false;
  }

  const hid_t attribute = H5Acreate2(id, name.c_str(), type, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  if (attribute < 0)
  {
    H5Sclose(dataspace);
    H5Tclose(type);
    return false;
  }

  const bool ok = H5Awrite(attribute, type, value.c_str()) >= 0;
  H5Aclose(attribute);
  H5Sclose(dataspace);
  H5Tclose(type);
  return ok;
}

bool
WriteDataset(hid_t id, const std::string& name, const Dataset& dataset)
{
  const hid_t dataspace =
    H5Screate_simple(static_cast<int>(dataset.dims.size()), dataset.dims.data(), nullptr);
  if (dataspace < 0)
    return false;

  const hid_t handle = H5Dcreate2(
    id, name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (handle < 0)
  {
    H5Sclose(dataspace);
    return false;
  }

  const bool ok =
    H5Dwrite(handle, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset.values.data()) >= 0;
  H5Dclose(handle);
  H5Sclose(dataspace);
  return ok;
}

/**
 * The datasets of the fissionable two-group set "set1" with two delayed groups, in the
 * prompt/delayed form. The delayed groups' shares of the delayed nu-fission depend on the
 * fissioning group: 1/2 and 1/2 in group 0, 1/3 and 2/3 in group 1.
 */
Datasets
PromptDelayedForm()
{
  return {
    {"inverse-velocity", Dataset1D({1.0, 0.5})},
    {"total", Dataset1D({1.0, 2.0})},
    {"fission", Dataset1D({0.5, 0.25})},
    {"nu-fission", Dataset1D({1.4, 1.0})},
    {"prompt-nu-fission", Dataset1D({1.0, 0.7})},
    {"chi", Dataset1D({0.6, 0.4})},
    {"chi-prompt", Dataset1D({0.7, 0.3})},
    {"decay-rate", Dataset1D({0.08, 0.2})},
    {"delayed-nu-fission", Dataset2D({0.2, 0.1, 0.2, 0.2}, 2, 2)},
    {"chi-delayed", Dataset2D({0.8, 0.2, 0.1, 0.9}, 2, 2)},
  };
}

/// The datasets of "set1" in the beta form, with delayed fractions independent of the group.
Datasets
BetaForm()
{
  return {
    {"inverse-velocity", Dataset1D({1.0, 0.5})},
    {"total", Dataset1D({1.0, 2.0})},
    {"fission", Dataset1D({0.5, 0.25})},
    {"nu-fission", Dataset1D({1.4, 1.0})},
    {"chi", Dataset1D({0.6, 0.4})},
    {"decay-rate", Dataset1D({0.08, 0.2})},
    {"beta", Dataset1D({0.002, 0.004})},
  };
}

/// Replaces, adds (a value) or removes (std::nullopt) datasets of `base`.
Datasets
With(Datasets base, const Datasets& changes)
{
  for (const auto& [name, dataset] : changes)
    base[name] = dataset;
  return base;
}

/**
 * Writes an OpenMC MGXS library with the two-group fissionable set "set1" made of `set1` and a
 * non-fissionable set "water", and returns its path. The file is rank-specific.
 */
std::string
WriteOpenMCFile(const Datasets& set1, const std::int64_t delayed_groups = 2)
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_openmc_delayed_test_r" + std::to_string(mpi_comm.rank()) + "_" +
                     std::to_string(now) + ".h5");

  const H5FileHandle file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  EXPECT_GE(file.Id(), 0);

  EXPECT_TRUE(WriteStringAttribute(file.Id(), "filetype", "mgxs"));
  EXPECT_TRUE(H5CreateAttribute<std::int64_t>(file.Id(), "energy_groups", 2));
  EXPECT_TRUE(H5CreateAttribute<std::int64_t>(file.Id(), "delayed_groups", delayed_groups));
  EXPECT_TRUE(H5WriteDataset1D<double>(file.Id(), "/group structure", {2.0e7, 1.0, 1.0e-5}));

  for (const std::string name : {"set1", "water"})
  {
    EXPECT_TRUE(H5CreateGroup(file.Id(), "/" + name));
    const hid_t group = H5Gopen2(file.Id(), ("/" + name).c_str(), H5P_DEFAULT);
    EXPECT_GE(group, 0);
    EXPECT_TRUE(H5CreateAttribute<unsigned int>(group, "order", 0));
    EXPECT_TRUE(H5CreateAttribute<bool>(group, "fissionable", name == "set1"));
    EXPECT_TRUE(WriteStringAttribute(group, "scatter_shape", "[G][G'][Order]"));
    H5Gclose(group);
  }

  EXPECT_TRUE(H5CreateGroup(file.Id(), "/set1/294K"));
  const hid_t xs = H5Gopen2(file.Id(), "/set1/294K", H5P_DEFAULT);
  EXPECT_GE(xs, 0);
  for (const auto& [name, dataset] : set1)
    if (dataset)
      EXPECT_TRUE(WriteDataset(xs, name, *dataset)) << name;
  H5Gclose(xs);

  EXPECT_TRUE(H5CreateGroup(file.Id(), "/water/294K"));
  const hid_t water_xs = H5Gopen2(file.Id(), "/water/294K", H5P_DEFAULT);
  EXPECT_GE(water_xs, 0);
  EXPECT_TRUE(H5WriteDataset1D<double>(water_xs, "inverse-velocity", {1.0, 0.5}));
  EXPECT_TRUE(H5WriteDataset1D<double>(water_xs, "total", {0.3, 0.4}));
  H5Gclose(water_xs);

  return path.string();
}

MultiGroupXS
Load(const Datasets& set1, const std::int64_t delayed_groups = 2)
{
  const auto file_name = WriteOpenMCFile(set1, delayed_groups);
  try
  {
    auto xs = MultiGroupXS::LoadFromOpenMC(file_name, "set1", 294.0);
    std::filesystem::remove(file_name);
    return xs;
  }
  catch (...)
  {
    std::filesystem::remove(file_name);
    throw;
  }
}

void
ExpectVectorNear(const std::vector<double>& actual, const std::vector<double>& expected)
{
  ASSERT_EQ(actual.size(), expected.size());
  for (size_t i = 0; i < expected.size(); ++i)
    EXPECT_NEAR(actual[i], expected[i], kTolerance) << "entry " << i;
}

/**
 * Checks the imported prompt production matrix chi_p (nu_p sigma_f)^T and the delayed production
 * sum_j chi_j gamma_j (nu_d sigma_f)^T against `expected_delayed`, [to group][from group].
 */
void
ExpectProduction(const MultiGroupXS& xs,
                 const std::vector<double>& chi_prompt,
                 const std::vector<double>& nu_prompt_sigma_f,
                 const std::vector<std::vector<double>>& expected_delayed)
{
  const auto& production = xs.GetProductionMatrix();
  ASSERT_EQ(production.size(), 2U);
  for (unsigned int g = 0; g < 2; ++g)
    for (unsigned int gp = 0; gp < 2; ++gp)
    {
      EXPECT_NEAR(production[g][gp], chi_prompt[g] * nu_prompt_sigma_f[gp], kTolerance);
      double delayed = 0.0;
      for (const auto& precursor : xs.GetPrecursors())
        delayed +=
          precursor.emission_spectrum[g] * precursor.fractional_yield * xs.GetNuDelayedSigmaF()[gp];
      EXPECT_NEAR(delayed, expected_delayed[g][gp], kTolerance) << g << " " << gp;
    }
}

} // namespace

TEST(OpenMCImportTest, ReadDelayedNeutronData)
{
  const auto xs = Load(PromptDelayedForm());

  EXPECT_TRUE(xs.IsFissionable());
  EXPECT_EQ(xs.GetNumGroups(), 2U);
  EXPECT_EQ(xs.GetNumPrecursors(), 2U);
  ExpectVectorNear(xs.GetNuSigmaF(), {1.4, 1.0});
  ExpectVectorNear(xs.GetNuPromptSigmaF(), {1.0, 0.7});
  ExpectVectorNear(xs.GetNuDelayedSigmaF(), {0.4, 0.3});

  // The yields are the shares of the delayed nu-fission summed over the groups: 0.3 and 0.4 of
  // 0.7. The delayed production is then sum_j chi_j gamma_j (nu_d sigma_f)^T.
  const auto& precursors = xs.GetPrecursors();
  ASSERT_EQ(precursors.size(), 2U);
  EXPECT_NEAR(precursors[0].decay_constant, 0.08, kTolerance);
  EXPECT_NEAR(precursors[1].decay_constant, 0.2, kTolerance);
  EXPECT_NEAR(precursors[0].fractional_yield, 0.3 / 0.7, kTolerance);
  EXPECT_NEAR(precursors[1].fractional_yield, 0.4 / 0.7, kTolerance);
  ExpectVectorNear(precursors[0].emission_spectrum, {0.8, 0.2});
  ExpectVectorNear(precursors[1].emission_spectrum, {0.1, 0.9});

  const double chi_d0 = (0.3 * 0.8 + 0.4 * 0.1) / 0.7;
  const double chi_d1 = (0.3 * 0.2 + 0.4 * 0.9) / 0.7;
  ExpectProduction(
    xs, {0.7, 0.3}, {1.0, 0.7}, {{chi_d0 * 0.4, chi_d0 * 0.3}, {chi_d1 * 0.4, chi_d1 * 0.3}});
}

TEST(OpenMCImportTest, DelayedProductionIsExactForGroupIndependentShares)
{
  // Shares 1/4 and 3/4 in both groups: sum_j chi_j (nu_d,j sigma_f)^T is reproduced exactly.
  const auto xs = Load(
    With(PromptDelayedForm(), {{"delayed-nu-fission", Dataset2D({0.1, 0.075, 0.3, 0.225}, 2, 2)}}));
  ASSERT_EQ(xs.GetPrecursors().size(), 2U);
  EXPECT_NEAR(xs.GetPrecursors()[0].fractional_yield, 0.25, kTolerance);
  EXPECT_NEAR(xs.GetPrecursors()[1].fractional_yield, 0.75, kTolerance);
  ExpectProduction(xs,
                   {0.7, 0.3},
                   {1.0, 0.7},
                   {{0.8 * 0.1 + 0.1 * 0.3, 0.8 * 0.075 + 0.1 * 0.225},
                    {0.2 * 0.1 + 0.9 * 0.3, 0.2 * 0.075 + 0.9 * 0.225}});
}

TEST(OpenMCImportTest, NuFissionIsPromptPlusDelayed)
{
  // As in OpenMC, the total is the sum of the parts, not the file's nu-fission.
  const auto xs = Load(With(PromptDelayedForm(), {{"nu-fission", Dataset1D({1.5, 0.9})}}));
  ExpectVectorNear(xs.GetNuSigmaF(), {1.4, 1.0});
}

TEST(OpenMCImportTest, PromptDelayedFormNeedsNeitherChiNorNuFission)
{
  const auto xs = Load(With(PromptDelayedForm(), {{"chi", std::nullopt}, {"nu-fission", {}}}));
  EXPECT_TRUE(xs.IsFissionable());
  EXPECT_TRUE(xs.GetChi().empty());
  ExpectVectorNear(xs.GetNuSigmaF(), {1.4, 1.0});
  ExpectVectorNear(xs.GetNuPromptSigmaF(), {1.0, 0.7});
}

TEST(OpenMCImportTest, PromptPartIsNuFissionMinusDelayedWithoutPromptNuFission)
{
  const auto xs = Load(With(PromptDelayedForm(), {{"prompt-nu-fission", std::nullopt}}));
  ExpectVectorNear(xs.GetNuSigmaF(), {1.4, 1.0});
  ExpectVectorNear(xs.GetNuPromptSigmaF(), {1.0, 0.7});
}

TEST(OpenMCImportTest, ChiIsPromptSpectrumWithoutChiPrompt)
{
  const auto xs = Load(With(PromptDelayedForm(), {{"chi-prompt", std::nullopt}}));
  const auto& production = xs.GetProductionMatrix();
  EXPECT_NEAR(production[0][0], 0.6 * 1.0, kTolerance);
  EXPECT_NEAR(production[1][1], 0.4 * 0.7, kTolerance);
}

TEST(OpenMCImportTest, SpectraAreNormalized)
{
  const auto xs = Load(With(PromptDelayedForm(),
                            {{"chi-prompt", Dataset1D({1.4, 0.6})},
                             {"chi-delayed", Dataset2D({4.0, 1.0, 0.5, 4.5}, 2, 2)}}));
  ExpectVectorNear(xs.GetPrecursors()[0].emission_spectrum, {0.8, 0.2});
  ExpectVectorNear(xs.GetPrecursors()[1].emission_spectrum, {0.1, 0.9});
  EXPECT_NEAR(xs.GetProductionMatrix()[0][0], 0.7 * 1.0, kTolerance);
  EXPECT_NEAR(xs.GetProductionMatrix()[1][0], 0.3 * 1.0, kTolerance);
}

TEST(OpenMCImportTest, OneDelayedGroupAcceptsVectors)
{
  const auto xs = Load(With(PromptDelayedForm(),
                            {{"delayed-nu-fission", Dataset1D({0.4, 0.3})},
                             {"chi-delayed", Dataset1D({0.5, 0.5})},
                             {"decay-rate", Dataset1D({0.1})}}),
                       1);
  ASSERT_EQ(xs.GetPrecursors().size(), 1U);
  EXPECT_NEAR(xs.GetPrecursors()[0].fractional_yield, 1.0, kTolerance);
  ExpectVectorNear(xs.GetNuDelayedSigmaF(), {0.4, 0.3});
}

TEST(OpenMCImportTest, BetaFormByDelayedGroup)
{
  const auto xs = Load(BetaForm());
  ASSERT_EQ(xs.GetPrecursors().size(), 2U);
  // delayed nu-fission beta_j nu_sigma_f, prompt (1 - 0.006) nu_sigma_f, all born with chi
  ExpectVectorNear(xs.GetNuSigmaF(), {1.4, 1.0});
  ExpectVectorNear(xs.GetNuDelayedSigmaF(), {0.006 * 1.4, 0.006 * 1.0});
  ExpectVectorNear(xs.GetNuPromptSigmaF(), {0.994 * 1.4, 0.994 * 1.0});
  EXPECT_NEAR(xs.GetPrecursors()[0].fractional_yield, 1.0 / 3.0, kTolerance);
  EXPECT_NEAR(xs.GetPrecursors()[1].fractional_yield, 2.0 / 3.0, kTolerance);
  for (const auto& precursor : xs.GetPrecursors())
    ExpectVectorNear(precursor.emission_spectrum, {0.6, 0.4});
  ExpectProduction(xs,
                   {0.6, 0.4},
                   {0.994 * 1.4, 0.994 * 1.0},
                   {{0.6 * 0.006 * 1.4, 0.6 * 0.006 * 1.0}, {0.4 * 0.006 * 1.4, 0.4 * 0.006}});
}

TEST(OpenMCImportTest, BetaFormByDelayedGroupAndGroup)
{
  // beta[j][g]: delayed fractions 0.003 + 0.003 in group 0 and 0.001 + 0.004 in group 1
  const auto xs = Load(With(BetaForm(), {{"beta", Dataset2D({0.003, 0.001, 0.003, 0.004}, 2, 2)}}));
  ExpectVectorNear(xs.GetNuDelayedSigmaF(), {0.006 * 1.4, 0.005 * 1.0});
  ExpectVectorNear(xs.GetNuPromptSigmaF(), {0.994 * 1.4, 0.995 * 1.0});
  const double total = 0.006 * 1.4 + 0.005;
  EXPECT_NEAR(xs.GetPrecursors()[0].fractional_yield, (0.003 * 1.4 + 0.001) / total, kTolerance);
  EXPECT_NEAR(xs.GetPrecursors()[1].fractional_yield, (0.003 * 1.4 + 0.004) / total, kTolerance);
  // All delayed neutrons are born with chi, so the delayed production is exact.
  ExpectProduction(xs,
                   {0.6, 0.4},
                   {0.994 * 1.4, 0.995 * 1.0},
                   {{0.6 * 0.006 * 1.4, 0.6 * 0.005}, {0.4 * 0.006 * 1.4, 0.4 * 0.005}});
}

TEST(OpenMCImportTest, BetaFormTakesPrecedence)
{
  // As in OpenMC, "beta" selects the beta form even when the prompt/delayed datasets exist.
  const auto xs = Load(With(PromptDelayedForm(), {{"beta", Dataset1D({0.002, 0.004})}}));
  ExpectVectorNear(xs.GetNuDelayedSigmaF(), {0.006 * 1.4, 0.006 * 1.0});
  for (const auto& precursor : xs.GetPrecursors())
    ExpectVectorNear(precursor.emission_spectrum, {0.6, 0.4});
}

TEST(OpenMCImportTest, NonFissionableMaterialHasNoPrecursors)
{
  const auto file_name = WriteOpenMCFile(PromptDelayedForm());
  const MultiGroupXS xs = MultiGroupXS::LoadFromOpenMC(file_name, "water", 294.0);
  std::filesystem::remove(file_name);

  EXPECT_FALSE(xs.IsFissionable());
  EXPECT_EQ(xs.GetNumGroups(), 2U);
  EXPECT_EQ(xs.GetNumPrecursors(), 0U);
  EXPECT_TRUE(xs.GetSigmaFission().empty());
  EXPECT_TRUE(xs.GetNuSigmaF().empty());
  EXPECT_TRUE(xs.GetNuPromptSigmaF().empty());
  EXPECT_TRUE(xs.GetNuDelayedSigmaF().empty());
  EXPECT_TRUE(xs.GetPrecursors().empty());
}

TEST(OpenMCImportTest, RejectsMalformedDelayedData)
{
  const std::map<std::string, Datasets> prompt_delayed_cases = {
    {"delayed-nu-fission with three delayed groups",
     {{"delayed-nu-fission", Dataset2D({0.2, 0.1, 0.2, 0.2, 0.0, 0.0}, 3, 2)}}},
    {"delayed-nu-fission with three groups",
     {{"delayed-nu-fission", Dataset2D({0.2, 0.1, 0.0, 0.2, 0.2, 0.0}, 2, 3)}}},
    {"one-dimensional delayed-nu-fission with two delayed groups",
     {{"delayed-nu-fission", Dataset1D({0.4, 0.3})}}},
    {"matrix delayed-nu-fission",
     {{"delayed-nu-fission", Dataset{{0.1, 0.1, 0.1, 0.1}, {1, 2, 2}}}}},
    {"negative delayed-nu-fission",
     {{"delayed-nu-fission", Dataset2D({0.2, -0.1, 0.2, 0.2}, 2, 2)}}},
    {"missing delayed-nu-fission", {{"delayed-nu-fission", std::nullopt}}},
    {"delayed above nu-fission",
     {{"prompt-nu-fission", std::nullopt}, {"nu-fission", Dataset1D({0.3, 1.0})}}},
    {"prompt-nu-fission with three groups", {{"prompt-nu-fission", Dataset1D({1.0, 0.7, 0.1})}}},
    {"chi-delayed with one delayed group", {{"chi-delayed", Dataset2D({0.8, 0.2}, 1, 2)}}},
    {"zero chi-delayed", {{"chi-delayed", Dataset2D({0.0, 0.0, 0.1, 0.9}, 2, 2)}}},
    {"missing chi-delayed", {{"chi-delayed", std::nullopt}}},
    {"zero chi-prompt", {{"chi-prompt", Dataset1D({0.0, 0.0})}}},
    {"no prompt spectrum", {{"chi-prompt", std::nullopt}, {"chi", std::nullopt}}},
    {"decay-rate with one value", {{"decay-rate", Dataset1D({0.08})}}},
    {"zero decay-rate", {{"decay-rate", Dataset1D({0.08, 0.0})}}},
    {"missing decay-rate", {{"decay-rate", std::nullopt}}},
  };
  for (const auto& [description, changes] : prompt_delayed_cases)
    EXPECT_THROW(Load(With(PromptDelayedForm(), changes)), std::logic_error) << description;

  const double inf = std::numeric_limits<double>::infinity();
  const std::map<std::string, Datasets> fission_cases = {
    {"fission with one group", {{"fission", Dataset1D({0.5})}}},
    {"infinite fission", {{"fission", Dataset1D({inf, 0.25})}}},
  };
  for (const auto& [description, changes] : fission_cases)
    EXPECT_THROW(Load(With(PromptDelayedForm(), changes)), std::logic_error) << description;

  const std::map<std::string, Datasets> beta_cases = {
    {"nu-fission with one group", {{"nu-fission", Dataset1D({1.4})}}},
    {"infinite nu-fission", {{"nu-fission", Dataset1D({inf, 1.0})}}},
    {"beta with one delayed group", {{"beta", Dataset1D({0.002})}}},
    {"beta summing above one", {{"beta", Dataset1D({0.6, 0.5})}}},
    {"negative beta", {{"beta", Dataset1D({0.002, -0.001})}}},
    {"beta with three groups", {{"beta", Dataset2D({0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, 2, 3)}}},
    {"beta without chi", {{"chi", std::nullopt}}},
  };
  for (const auto& [description, changes] : beta_cases)
    EXPECT_THROW(Load(With(BetaForm(), changes)), std::logic_error) << description;

  EXPECT_THROW(Load(PromptDelayedForm(), -1), std::logic_error);
}

TEST(OpenMCImportTest, DelayedGroupsAttribute)
{
  // Without the attribute the library has no delayed groups: nu-fission is all prompt.
  auto file_name = WriteOpenMCFile(PromptDelayedForm());
  {
    const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    ASSERT_GE(H5Adelete(file.Id(), "delayed_groups"), 0);
  }
  const auto xs = MultiGroupXS::LoadFromOpenMC(file_name, "set1", 294.0);
  EXPECT_EQ(xs.GetNumPrecursors(), 0U);
  ExpectVectorNear(xs.GetNuSigmaF(), {1.4, 1.0});
  EXPECT_NEAR(xs.GetProductionMatrix()[0][0], 0.6 * 1.4, kTolerance);

  // An attribute that cannot be read as an integer is rejected rather than taken as zero.
  {
    const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    ASSERT_TRUE(WriteStringAttribute(file.Id(), "delayed_groups", "two"));
  }
  EXPECT_THROW(MultiGroupXS::LoadFromOpenMC(file_name, "set1", 294.0), std::logic_error);
  std::filesystem::remove(file_name);
}
