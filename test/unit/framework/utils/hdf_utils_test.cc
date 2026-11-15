// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/runtime.h"
#include "framework/utils/hdf_utils.h"
#include "gtest/gtest.h"
#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
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
