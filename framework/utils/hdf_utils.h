// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "hdf5.h"
#include <vector>
#include <string>
#include <type_traits>

namespace opensn
{
template <typename T>
hid_t get_datatype();

template <>
inline hid_t
get_datatype<char>()
{
  return H5T_NATIVE_CHAR;
}

template <>
inline hid_t
get_datatype<signed char>()
{
  return H5T_NATIVE_SCHAR;
}

template <>
inline hid_t
get_datatype<unsigned char>()
{
  return H5T_NATIVE_UCHAR;
}

template <>
inline hid_t
get_datatype<short>()
{
  return H5T_NATIVE_SHORT;
}

template <>
inline hid_t
get_datatype<unsigned short>()
{
  return H5T_NATIVE_USHORT;
}

template <>
inline hid_t
get_datatype<int>()
{
  return H5T_NATIVE_INT;
}

template <>
inline hid_t
get_datatype<unsigned int>()
{
  return H5T_NATIVE_UINT;
}

template <>
inline hid_t
get_datatype<long>()
{
  return H5T_NATIVE_LONG;
}

template <>
inline hid_t
get_datatype<unsigned long>()
{
  return H5T_NATIVE_ULONG;
}

template <>
inline hid_t
get_datatype<long long>()
{
  return H5T_NATIVE_LLONG;
}

template <>
inline hid_t
get_datatype<unsigned long long>()
{
  return H5T_NATIVE_ULLONG;
}

template <>
inline hid_t
get_datatype<float>()
{
  return H5T_NATIVE_FLOAT;
}

template <>
inline hid_t
get_datatype<double>()
{
  return H5T_NATIVE_DOUBLE;
}

template <>
inline hid_t
get_datatype<bool>()
{
  return H5T_NATIVE_UCHAR;
}

inline bool
H5Has(hid_t id, const std::string& name)
{
  return (H5Lexists(id, name.c_str(), H5P_DEFAULT) > 0);
}

inline bool
H5CreateGroup(hid_t id, const std::string& name)
{
  hid_t group = H5Gcreate2(id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (group < 0)
    return false;
  H5Gclose(group);
  return true;
}

template <typename T>
bool
H5WriteDataset1D(hid_t id, const std::string& name, const std::vector<T>& data)
{
  static_assert(
    !std::is_same_v<T, bool>,
    "H5WriteDataset1D does not support std::vector<bool>. Use std::vector<unsigned char>.");

  bool success = false;

  hsize_t dim[1] = {static_cast<hsize_t>(data.size())};
  auto dataspace = H5Screate_simple(1, dim, NULL);
  if (dataspace != H5I_INVALID_HID)
  {
    auto dataset = H5Dcreate2(
      id, name.c_str(), get_datatype<T>(), dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset != H5I_INVALID_HID)
    {
      if (data.empty())
        success = true;
      else
      {
        if (H5Dwrite(dataset, get_datatype<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) >= 0)
          success = true;
      }
      H5Dclose(dataset);
    }
    H5Sclose(dataspace);
  }

  return success;
}

template <typename T>
bool
H5CreateAttribute(hid_t id, const std::string& name, const T& data)
{
  static_assert(!std::is_same_v<T, std::string>, "H5CreateAttribute does not support std::string.");

  bool success = false;

  auto dataspace = H5Screate(H5S_SCALAR);
  if (dataspace != H5I_INVALID_HID)
  {
    auto attribute =
      H5Acreate2(id, name.c_str(), get_datatype<T>(), dataspace, H5P_DEFAULT, H5P_DEFAULT);
    if (attribute != H5I_INVALID_HID)
    {
      if (H5Awrite(attribute, get_datatype<T>(), &data) >= 0)
        success = true;
      H5Aclose(attribute);
    }
    H5Sclose(dataspace);
  }

  return success;
}

template <typename T>
bool
H5ReadDataset1D(hid_t id, const std::string& name, std::vector<T>& data)
{
  static_assert(
    !std::is_same_v<T, bool>,
    "H5ReadDataset1D does not support std::vector<bool>. Use std::vector<unsigned char>.");

  bool success = false;

  data.clear();

  auto dataset = H5Dopen2(id, name.c_str(), H5P_DEFAULT);
  if (dataset != H5I_INVALID_HID)
  {
    auto dataspace = H5Dget_space(dataset);
    if (dataspace != H5I_INVALID_HID)
    {
      const int rank = H5Sget_simple_extent_ndims(dataspace);
      if (rank == 1)
      {
        const auto npoints = H5Sget_simple_extent_npoints(dataspace);
        if (npoints >= 0)
        {
          data.resize(static_cast<size_t>(npoints));
          if (npoints == 0 or
              H5Dread(dataset, get_datatype<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) >= 0)
          {
            success = true;
          }
          else
          {
            data.clear();
            data.shrink_to_fit();
          }
        }
      }
      H5Sclose(dataspace);
    }
    H5Dclose(dataset);

    if (success)
      return true;
  }

  if (H5Aexists(id, name.c_str()) > 0)
  {
    auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (attribute != H5I_INVALID_HID)
    {
      auto aspace = H5Aget_space(attribute);
      if (aspace != H5I_INVALID_HID)
      {
        const auto npoints = H5Sget_simple_extent_npoints(aspace);
        if (npoints >= 0)
        {
          data.resize(static_cast<size_t>(npoints));
          if (npoints == 0 or H5Aread(attribute, get_datatype<T>(), data.data()) >= 0)
          {
            success = true;
          }
          else
          {
            data.clear();
            data.shrink_to_fit();
          }
        }
        H5Sclose(aspace);
      }
      H5Aclose(attribute);
    }
  }

  return success;
}

template <typename T>
bool
H5ReadAttribute(hid_t id, const std::string& name, T& value)
{
  bool success = false;

  if (H5Aexists(id, name.c_str()) > 0)
  {
    auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (attribute != H5I_INVALID_HID)
    {
      if (H5Aread(attribute, get_datatype<T>(), &value) >= 0)
        success = true;
      H5Aclose(attribute);
    }
  }

  return success;
}

template <>
inline bool
H5ReadAttribute<std::string>(hid_t id, const std::string& name, std::string& value)
{
  // NOTE: This implementation assumes the attribute was written as a
  // fixed-length C string and pads with a NUL terminator.

  bool success = false;

  if (H5Aexists(id, name.c_str()) <= 0)
    return false;

  auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
  if (attribute == H5I_INVALID_HID)
    return false;

  auto size = static_cast<size_t>(H5Aget_storage_size(attribute));
  hid_t string_type = H5Tcopy(H5T_C_S1);
  if (string_type >= 0)
  {
    if (size == 0)
    {
      value.clear();
      success = true;
    }
    else
    {
      H5Tset_size(string_type, size + 1);
      H5Tset_strpad(string_type, H5T_STR_NULLTERM);
      std::vector<char> buffer(size + 1, '\0');
      if (H5Aread(attribute, string_type, buffer.data()) >= 0)
      {
        value = buffer.data();
        success = true;
      }
    }
    H5Tclose(string_type);
  }
  H5Aclose(attribute);

  return success;
}

template <>
inline bool
H5ReadAttribute<bool>(hid_t id, const std::string& name, bool& value)
{
  bool success = false;

  if (H5Aexists(id, name.c_str()) <= 0)
    return false;

  auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
  if (attribute == H5I_INVALID_HID)
    return false;

  unsigned char tmp = 0;
  if (H5Aread(attribute, H5T_NATIVE_UCHAR, &tmp) >= 0)
  {
    value = static_cast<bool>(tmp);
    success = true;
  }
  H5Aclose(attribute);

  return success;
}

template <typename T>
bool
H5ReadGroupAttribute(hid_t id, const std::string& group_id, const std::string& name, T& value)
{
  bool success = false;

  auto group = H5Gopen2(id, group_id.c_str(), H5P_DEFAULT);
  if (group != H5I_INVALID_HID)
  {
    success = H5ReadAttribute<T>(group, name, value);
    H5Gclose(group);
  }

  return success;
}

template <>
inline bool
H5ReadGroupAttribute<std::string>(hid_t id,
                                  const std::string& group_id,
                                  const std::string& name,
                                  std::string& value)
{
  bool success = false;

  auto group = H5Gopen2(id, group_id.c_str(), H5P_DEFAULT);
  if (group != H5I_INVALID_HID)
  {
    success = H5ReadAttribute<std::string>(group, name, value);
    H5Gclose(group);
  }

  return success;
}

template <>
inline bool
H5ReadGroupAttribute<bool>(hid_t id,
                           const std::string& group_id,
                           const std::string& name,
                           bool& value)
{
  bool success = false;

  auto group = H5Gopen2(id, group_id.c_str(), H5P_DEFAULT);
  if (group != H5I_INVALID_HID)
  {
    success = H5ReadAttribute<bool>(group, name, value);
    H5Gclose(group);
  }

  return success;
}

} // namespace opensn
