// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "hdf5.h"
#include <stdexcept>
#include <vector>
#include <string>
#include <cassert>

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

inline bool
H5Has(hid_t id, const std::string& name)
{
  return (H5Lexists(id, name.c_str(), H5P_DEFAULT) > 0);
}

inline bool
H5CreateGroup(hid_t id, const std::string& name)
{
  return (H5Gcreate2(id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) >= 0);
}

template <typename T>
bool
H5WriteDataset1D(hid_t id, const std::string& name, std::vector<T> data)
{
  bool success = false;

  hsize_t dim[1], maxdims[1];
  dim[0] = data.size();
  maxdims[0] = data.size();
  auto dataspace = H5Screate_simple(1, dim, maxdims);
  if (dataspace != H5I_INVALID_HID)
  {
    auto dataset = H5Dcreate2(
      id, name.c_str(), get_datatype<T>(), dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset != H5I_INVALID_HID)
    {
      if (H5Dwrite(dataset, get_datatype<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) >= 0)
        success = true;
      H5Dclose(dataset);
    }
    H5Sclose(dataspace);
  }

  return success;
}

template <typename T>
bool
H5CreateAttribute(hid_t id, const std::string& name, T& data)
{
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
  bool success = true;

  auto dataset = H5Dopen2(id, name.c_str(), H5P_DEFAULT);
  if (dataset != H5I_INVALID_HID)
  {
    auto dataspace = H5Dget_space(dataset);
    if (dataspace != H5I_INVALID_HID)
    {
      hsize_t dims[1];
      if (H5Sget_simple_extent_dims(dataspace, dims, NULL) == 1)
      {
        data.resize(dims[0]);
        if (H5Dread(dataset, get_datatype<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0)
        {
          data.clear();
          data.shrink_to_fit();
          success = false;
        }
      }
      H5Sclose(dataspace);
    }
    H5Dclose(dataset);
  }
  else
  {
    if (H5Aexists(id, name.c_str()))
    {
      auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
      if (attribute != H5I_INVALID_HID)
      {
        size_t size = static_cast<size_t>(H5Aget_storage_size(attribute));
        if (size > 0)
        {
          size_t num_elements = size / sizeof(T);
          data.resize(num_elements);
          if (H5Aread(attribute, get_datatype<T>(), data.data()) < 0)
          {
            data.clear();
            data.shrink_to_fit();
            success = false;
          }
        }
        H5Aclose(attribute);
      }
    }
  }

  return success;
}

template <typename T>
bool
H5ReadAttribute(hid_t id, const std::string& name, T& value)
{
  bool success = false;

  if (H5Aexists(id, name.c_str()))
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
  bool success = false;

  if (H5Aexists(id, name.c_str()))
  {
    auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (attribute != H5I_INVALID_HID)
    {
      size_t size = static_cast<size_t>(H5Aget_storage_size(attribute));
      if (size > 0)
      {
        hid_t string_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(string_type, size + 1);
        std::vector<char> buffer(size + 1);
        if (H5Aread(attribute, string_type, buffer.data()) >= 0)
        {
          value = buffer.data();
          success = true;
        }
      }
      H5Aclose(attribute);
    }
  }

  return success;
}

template <>
inline bool
H5ReadAttribute<bool>(hid_t id, const std::string& name, bool& value)
{
  bool success = false;

  if (H5Aexists(id, name.c_str()))
  {
    auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (attribute != H5I_INVALID_HID)
    {
      if (H5Aread(attribute, H5T_NATIVE_INT, &value) >= 0)
        success = true;
      H5Aclose(attribute);
    }
  }

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
