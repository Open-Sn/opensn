#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/utils/hdf_utils.h"
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>
#include <gtest/gtest.h>

using namespace opensn;

namespace
{

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
WriteDataset2D(hid_t id,
               const std::string& name,
               const std::vector<double>& data,
               const hsize_t dim0,
               const hsize_t dim1)
{
  const hsize_t dims[2] = {dim0, dim1};
  const hid_t dataspace = H5Screate_simple(2, dims, nullptr);
  if (dataspace < 0)
    return false;

  const hid_t dataset = H5Dcreate2(
    id, name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0)
  {
    H5Sclose(dataspace);
    return false;
  }

  const bool ok =
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) >= 0;
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return ok;
}

std::string
MakeTempOpenMCDelayedXSFile()
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_openmc_delayed_test_" + std::to_string(now) + ".h5");

  const H5FileHandle file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  EXPECT_GE(file.Id(), 0);

  EXPECT_TRUE(WriteStringAttribute(file.Id(), "filetype", "mgxs"));
  EXPECT_TRUE(H5CreateAttribute<std::int64_t>(file.Id(), "energy_groups", 2));
  EXPECT_TRUE(H5CreateAttribute<std::int64_t>(file.Id(), "delayed_groups", 2));
  EXPECT_TRUE(H5WriteDataset1D<double>(file.Id(), "/group structure", {2.0e7, 1.0, 1.0e-5}));

  EXPECT_TRUE(H5CreateGroup(file.Id(), "/set1"));
  const hid_t set1 = H5Gopen2(file.Id(), "/set1", H5P_DEFAULT);
  EXPECT_GE(set1, 0);
  EXPECT_TRUE(H5CreateAttribute<unsigned int>(set1, "order", 0));
  EXPECT_TRUE(H5CreateAttribute<bool>(set1, "fissionable", true));
  EXPECT_TRUE(WriteStringAttribute(set1, "scatter_shape", "[G][G'][Order]"));
  H5Gclose(set1);

  EXPECT_TRUE(H5CreateGroup(file.Id(), "/set1/294K"));
  const hid_t xs = H5Gopen2(file.Id(), "/set1/294K", H5P_DEFAULT);
  EXPECT_GE(xs, 0);

  EXPECT_TRUE(H5WriteDataset1D<double>(xs, "inverse-velocity", {1.0, 0.5}));
  EXPECT_TRUE(H5WriteDataset1D<double>(xs, "total", {1.0, 2.0}));
  EXPECT_TRUE(H5WriteDataset1D<double>(xs, "fission", {0.5, 0.25}));
  EXPECT_TRUE(H5WriteDataset1D<double>(xs, "nu-fission", {1.4, 1.0}));
  EXPECT_TRUE(H5WriteDataset1D<double>(xs, "prompt-nu-fission", {1.0, 0.7}));
  EXPECT_TRUE(H5WriteDataset1D<double>(xs, "chi", {0.6, 0.4}));
  EXPECT_TRUE(H5WriteDataset1D<double>(xs, "chi-prompt", {0.7, 0.3}));
  EXPECT_TRUE(H5WriteDataset1D<double>(xs, "decay-rate", {0.08, 0.2}));
  EXPECT_TRUE(WriteDataset2D(xs, "delayed-nu-fission", {0.2, 0.1, 0.2, 0.2}, 2, 2));
  EXPECT_TRUE(WriteDataset2D(xs, "chi-delayed", {0.8, 0.2, 0.1, 0.9}, 2, 2));

  H5Gclose(xs);

  EXPECT_TRUE(H5CreateGroup(file.Id(), "/water"));
  const hid_t water = H5Gopen2(file.Id(), "/water", H5P_DEFAULT);
  EXPECT_GE(water, 0);
  EXPECT_TRUE(H5CreateAttribute<unsigned int>(water, "order", 0));
  EXPECT_TRUE(H5CreateAttribute<bool>(water, "fissionable", false));
  EXPECT_TRUE(WriteStringAttribute(water, "scatter_shape", "[G][G'][Order]"));
  H5Gclose(water);

  EXPECT_TRUE(H5CreateGroup(file.Id(), "/water/294K"));
  const hid_t water_xs = H5Gopen2(file.Id(), "/water/294K", H5P_DEFAULT);
  EXPECT_GE(water_xs, 0);
  EXPECT_TRUE(H5WriteDataset1D<double>(water_xs, "inverse-velocity", {1.0, 0.5}));
  EXPECT_TRUE(H5WriteDataset1D<double>(water_xs, "total", {0.3, 0.4}));
  H5Gclose(water_xs);

  return path.string();
}

} // namespace

TEST(OpenMCImportTest, ReadDelayedNeutronData)
{
  const std::string fname = MakeTempOpenMCDelayedXSFile();
  const MultiGroupXS xs = MultiGroupXS::LoadFromOpenMC(fname, "set1", 294.0);

  EXPECT_TRUE(xs.IsFissionable());
  EXPECT_EQ(xs.GetNumGroups(), 2u);
  EXPECT_EQ(xs.GetNumPrecursors(), 2u);

  ASSERT_EQ(xs.GetNuSigmaF().size(), 2u);
  EXPECT_NEAR(xs.GetNuSigmaF()[0], 1.4, 1.0e-12);
  EXPECT_NEAR(xs.GetNuSigmaF()[1], 1.0, 1.0e-12);

  ASSERT_EQ(xs.GetNuPromptSigmaF().size(), 2u);
  EXPECT_NEAR(xs.GetNuPromptSigmaF()[0], 1.0, 1.0e-12);
  EXPECT_NEAR(xs.GetNuPromptSigmaF()[1], 0.7, 1.0e-12);

  ASSERT_EQ(xs.GetNuDelayedSigmaF().size(), 2u);
  EXPECT_NEAR(xs.GetNuDelayedSigmaF()[0], 0.4, 1.0e-12);
  EXPECT_NEAR(xs.GetNuDelayedSigmaF()[1], 0.3, 1.0e-12);

  const auto& production = xs.GetProductionMatrix();
  ASSERT_EQ(production.size(), 2u);
  EXPECT_NEAR(production[0][0], 0.7, 1.0e-12);
  EXPECT_NEAR(production[1][0], 0.3, 1.0e-12);
  EXPECT_NEAR(production[0][1], 0.49, 1.0e-12);
  EXPECT_NEAR(production[1][1], 0.21, 1.0e-12);

  const auto& precursors = xs.GetPrecursors();
  ASSERT_EQ(precursors.size(), 2u);
  EXPECT_NEAR(precursors[0].decay_constant, 0.08, 1.0e-12);
  EXPECT_NEAR(precursors[1].decay_constant, 0.2, 1.0e-12);
  EXPECT_NEAR(precursors[0].fractional_yield, 0.3 / 0.7, 1.0e-12);
  EXPECT_NEAR(precursors[1].fractional_yield, 0.4 / 0.7, 1.0e-12);
  ASSERT_EQ(precursors[0].emission_spectrum.size(), 2u);
  EXPECT_NEAR(precursors[0].emission_spectrum[0], 0.8, 1.0e-12);
  EXPECT_NEAR(precursors[0].emission_spectrum[1], 0.2, 1.0e-12);
  EXPECT_NEAR(precursors[1].emission_spectrum[0], 0.1, 1.0e-12);
  EXPECT_NEAR(precursors[1].emission_spectrum[1], 0.9, 1.0e-12);

  std::filesystem::remove(fname);
}

TEST(OpenMCImportTest, NonFissionableMaterialHasNoPrecursors)
{
  const std::string fname = MakeTempOpenMCDelayedXSFile();
  const MultiGroupXS xs = MultiGroupXS::LoadFromOpenMC(fname, "water", 294.0);

  EXPECT_FALSE(xs.IsFissionable());
  EXPECT_EQ(xs.GetNumGroups(), 2u);
  EXPECT_EQ(xs.GetNumPrecursors(), 0u);
  EXPECT_TRUE(xs.GetSigmaFission().empty());
  EXPECT_TRUE(xs.GetNuSigmaF().empty());
  EXPECT_TRUE(xs.GetNuPromptSigmaF().empty());
  EXPECT_TRUE(xs.GetNuDelayedSigmaF().empty());
  EXPECT_TRUE(xs.GetPrecursors().empty());

  std::filesystem::remove(fname);
}
