// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/runtime.h"
#include "framework/utils/hdf_utils.h"
#include "gtest/gtest.h"
#include <array>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

using namespace opensn;

namespace
{

/// Write an attribute of the given number of unsigned int elements.
void
WriteArrayAttribute(hid_t id, const std::string& name, const std::array<unsigned int, 3>& values)
{
  const hsize_t size = values.size();
  const hid_t dataspace = H5Screate_simple(1, &size, nullptr);
  ASSERT_NE(dataspace, H5I_INVALID_HID);
  const hid_t attribute =
    H5Acreate2(id, name.c_str(), H5T_NATIVE_UINT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_NE(attribute, H5I_INVALID_HID);
  EXPECT_GE(H5Awrite(attribute, H5T_NATIVE_UINT, values.data()), 0);
  H5Aclose(attribute);
  H5Sclose(dataspace);
}

/// Write a one-element attribute with a simple, not scalar, dataspace.
void
WriteOneElementAttribute(hid_t id, const std::string& name, unsigned int value)
{
  const hsize_t size = 1;
  const hid_t dataspace = H5Screate_simple(1, &size, nullptr);
  ASSERT_NE(dataspace, H5I_INVALID_HID);
  const hid_t attribute =
    H5Acreate2(id, name.c_str(), H5T_NATIVE_UINT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_NE(attribute, H5I_INVALID_HID);
  EXPECT_GE(H5Awrite(attribute, H5T_NATIVE_UINT, &value), 0);
  H5Aclose(attribute);
  H5Sclose(dataspace);
}

/// Write an attribute of doubles with a simple dataspace.
void
WriteDoubleAttribute(hid_t id, const std::string& name, const std::vector<double>& values)
{
  const hsize_t size = values.size();
  const hid_t dataspace = H5Screate_simple(1, &size, nullptr);
  ASSERT_NE(dataspace, H5I_INVALID_HID);
  const hid_t attribute =
    H5Acreate2(id, name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  ASSERT_NE(attribute, H5I_INVALID_HID);
  EXPECT_GE(H5Awrite(attribute, H5T_NATIVE_DOUBLE, values.data()), 0);
  H5Aclose(attribute);
  H5Sclose(dataspace);
}

} // namespace

/// A scalar attribute read must reject an attribute with more than one element instead of writing
/// past the end of the value.
TEST(HDFUtils, ReadAttributeRejectsArrays)
{
  const auto file_name =
    (std::filesystem::temp_directory_path() / ("opensn_hdf_utils_" + std::to_string(getpid()) +
                                               "_r" + std::to_string(mpi_comm.rank()) + ".h5"))
      .string();
  {
    const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_GE(file.Id(), 0);
    ASSERT_TRUE(H5CreateAttribute<unsigned int>(file.Id(), "scalar", 7U));
    WriteOneElementAttribute(file.Id(), "one", 5U);
    WriteArrayAttribute(file.Id(), "array", {1U, 2U, 3U});
  }
  {
    const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    ASSERT_GE(file.Id(), 0);

    unsigned int value = 0;
    EXPECT_TRUE(H5ReadAttribute<unsigned int>(file.Id(), "scalar", value));
    EXPECT_EQ(value, 7U);
    EXPECT_TRUE(H5ReadAttribute<unsigned int>(file.Id(), "one", value));
    EXPECT_EQ(value, 5U);

    // The value is the first element of a buffer of three, so a read of the whole attribute stays
    // in bounds and is detected by the return value and the untouched elements.
    std::array<unsigned int, 3> buffer = {0U, 0U, 0U};
    EXPECT_FALSE(H5ReadAttribute<unsigned int>(file.Id(), "array", buffer[0]));
    EXPECT_EQ(buffer, (std::array<unsigned int, 3>{0U, 0U, 0U}));
    EXPECT_FALSE(H5ReadOptionalAttribute<unsigned int>(file.Id(), "array", buffer[0]));
    EXPECT_EQ(buffer, (std::array<unsigned int, 3>{0U, 0U, 0U}));

    bool flag = false;
    EXPECT_FALSE(H5ReadAttribute<bool>(file.Id(), "array", flag));

    // A missing optional attribute leaves the value unchanged.
    value = 11U;
    EXPECT_TRUE(H5ReadOptionalAttribute<unsigned int>(file.Id(), "missing", value));
    EXPECT_EQ(value, 11U);
  }
  std::filesystem::remove(file_name);
}

/// A dataset or attribute read with an expected size must reject any other extent, leave the
/// output empty, and close the datasets and attributes it opened. A dataset of the wrong size must
/// not be replaced by an attribute of the same name. (HDF5 1.10 cannot count open dataspaces, so
/// their release is not checked here.)
TEST(HDFUtils, ReadDataset1DChecksTheExpectedSize)
{
  const auto file_name =
    (std::filesystem::temp_directory_path() / ("opensn_hdf_utils_size_" + std::to_string(getpid()) +
                                               "_r" + std::to_string(mpi_comm.rank()) + ".h5"))
      .string();
  const std::vector<double> values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  {
    const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_GE(file.Id(), 0);
    ASSERT_TRUE(H5WriteDataset1D(file.Id(), "values", values));
    WriteArrayAttribute(file.Id(), "array", {1U, 2U, 3U});
    // A 7-value dataset and an 8-value attribute with the same name.
    ASSERT_TRUE(H5WriteDataset1D(file.Id(), "both", std::vector<double>(7, 1.0)));
    WriteDoubleAttribute(file.Id(), "both", std::vector<double>(8, 2.0));
  }
  {
    const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    ASSERT_GE(file.Id(), 0);

    std::vector<double> data = {-1.0};
    EXPECT_FALSE(H5ReadDataset1D(file.Id(), "values", data, std::uint64_t{3}));
    EXPECT_TRUE(data.empty());
    EXPECT_TRUE(H5ReadDataset1D(file.Id(), "values", data, std::uint64_t{10}));
    EXPECT_EQ(data, values);
    EXPECT_TRUE(H5ReadDataset1D(file.Id(), "values", data));
    EXPECT_EQ(data, values);

    std::vector<unsigned int> array = {7U};
    EXPECT_FALSE(H5ReadDataset1D(file.Id(), "array", array, std::uint64_t{2}));
    EXPECT_TRUE(array.empty());
    EXPECT_TRUE(H5ReadDataset1D(file.Id(), "array", array, std::uint64_t{3}));
    EXPECT_EQ(array, (std::vector<unsigned int>{1U, 2U, 3U}));

    std::vector<double> both = {-1.0};
    EXPECT_FALSE(H5ReadDataset1D(file.Id(), "both", both, std::uint64_t{8}));
    EXPECT_TRUE(both.empty());
    EXPECT_TRUE(H5ReadDataset1D(file.Id(), "both", both, std::uint64_t{7}));
    EXPECT_EQ(both, std::vector<double>(7, 1.0));

    // No dataset or attribute of the file is still open.
    EXPECT_EQ(H5Fget_obj_count(file.Id(), H5F_OBJ_DATASET | H5F_OBJ_LOCAL), 0);
    EXPECT_EQ(H5Fget_obj_count(file.Id(), H5F_OBJ_ATTR | H5F_OBJ_LOCAL), 0);
  }
  std::filesystem::remove(file_name);
}

/// A read whose buffer cannot be allocated throws, and closes the dataset it opened.
TEST(HDFUtils, ReadDataset1DClosesTheDatasetWhenAllocationFails)
{
  const auto file_name =
    (std::filesystem::temp_directory_path() / ("opensn_hdf_utils_huge_" + std::to_string(getpid()) +
                                               "_r" + std::to_string(mpi_comm.rank()) + ".h5"))
      .string();
  {
    const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_GE(file.Id(), 0);
    // A chunked dataset that is never written stores nothing, whatever its extent. Its 2^61
    // elements exceed what a vector of doubles can hold, so resizing throws before allocating.
    // The identifiers are guarded as soon as they are made, so a failed check releases them.
    const hsize_t size = hsize_t{1} << 61;
    const hsize_t chunk = 1024;
    {
      const H5IdGuard dataspace(H5Screate_simple(1, &size, nullptr), H5Sclose);
      const H5IdGuard properties(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
      ASSERT_GE(dataspace.Id(), 0);
      ASSERT_GE(properties.Id(), 0);
      ASSERT_GE(H5Pset_chunk(properties.Id(), 1, &chunk), 0);
      const H5IdGuard dataset(H5Dcreate2(file.Id(),
                                         "huge",
                                         H5T_NATIVE_UCHAR,
                                         dataspace.Id(),
                                         H5P_DEFAULT,
                                         properties.Id(),
                                         H5P_DEFAULT),
                              H5Dclose);
      ASSERT_GE(dataset.Id(), 0);
    }

    std::vector<double> data;
    EXPECT_THROW(H5ReadDataset1D<double>(file.Id(), "huge", data), std::length_error);
    EXPECT_EQ(H5Fget_obj_count(file.Id(), H5F_OBJ_DATASET | H5F_OBJ_LOCAL), 0);
  }
  std::filesystem::remove(file_name);
}
