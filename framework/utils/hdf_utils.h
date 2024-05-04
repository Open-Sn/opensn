// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "hdf5.h"
#include <vector>

namespace opensn
{
template <typename T>
hid_t get_datatype();

template <>
hid_t
get_datatype<char>()
{
  return H5T_NATIVE_CHAR;
}

template <>
hid_t
get_datatype<signed char>()
{
  return H5T_NATIVE_SCHAR;
}

template <>
hid_t
get_datatype<unsigned char>()
{
  return H5T_NATIVE_UCHAR;
}

template <>
hid_t
get_datatype<short>()
{
  return H5T_NATIVE_SHORT;
}

template <>
hid_t
get_datatype<unsigned short>()
{
  return H5T_NATIVE_USHORT;
}

template <>
hid_t
get_datatype<int>()
{
  return H5T_NATIVE_INT;
}

template <>
hid_t
get_datatype<unsigned int>()
{
  return H5T_NATIVE_UINT;
}

template <>
hid_t
get_datatype<long>()
{
  return H5T_NATIVE_LONG;
}

template <>
hid_t
get_datatype<unsigned long>()
{
  return H5T_NATIVE_ULONG;
}

template <>
hid_t
get_datatype<long long>()
{
  return H5T_NATIVE_LLONG;
}

template <>
hid_t
get_datatype<unsigned long long>()
{
  return H5T_NATIVE_ULLONG;
}

template <>
hid_t
get_datatype<float>()
{
  return H5T_NATIVE_FLOAT;
}

template <>
hid_t
get_datatype<double>()
{
  return H5T_NATIVE_DOUBLE;
}

bool
H5Has(hid_t id, const std::string& name)
{
  return (H5Lexists(id, name.c_str(), H5P_DEFAULT) > 0);
}

template <typename T>
std::vector<T>
H5ReadDataset1D(hid_t id, const std::string& name)
{
  std::vector<T> data;

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
          }
        }
        H5Aclose(attribute);
      }
    }
  }

  return data;
}

template <typename T>
bool
H5ReadAttribute(hid_t id, const std::string& name, T& value)
{
  bool retval = false;

  if (H5Aexists(id, name.c_str()))
  {
    auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (attribute != H5I_INVALID_HID)
    {
      if (H5Aread(attribute, get_datatype<T>(), &value) >= 0)
        retval = true;
      H5Aclose(attribute);
    }
  }

  return retval;
}

template <>
bool
H5ReadAttribute<std::string>(hid_t id, const std::string& name, std::string& value)
{
  bool retval = false;

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
          retval = true;
        }
      }
      H5Aclose(attribute);
    }
  }

  return retval;
}

template <>
bool
H5ReadAttribute<bool>(hid_t id, const std::string& name, bool& value)
{
  bool retval = false;

  if (H5Aexists(id, name.c_str()))
  {
    auto attribute = H5Aopen(id, name.c_str(), H5P_DEFAULT);
    if (attribute != H5I_INVALID_HID)
    {
      if (H5Aread(attribute, H5T_NATIVE_INT, &value) >= 0)
        retval = true;
      H5Aclose(attribute);
    }
  }

  return retval;
}

template <typename T>
bool
H5ReadGroupAttribute(hid_t id, const std::string& group_id, const std::string& name, T& value)
{
  bool retval = false;

  auto group = H5Gopen2(id, group_id.c_str(), H5P_DEFAULT);
  if (group != H5I_INVALID_HID)
  {
    retval = H5ReadAttribute<T>(group, name, value);
    H5Gclose(group);
  }

  return retval;
}

template <>
bool
H5ReadGroupAttribute<std::string>(hid_t id,
                                  const std::string& group_id,
                                  const std::string& name,
                                  std::string& value)
{
  bool retval = false;

  auto group = H5Gopen2(id, group_id.c_str(), H5P_DEFAULT);
  if (group != H5I_INVALID_HID)
  {
    retval = H5ReadAttribute<std::string>(group, name, value);
    H5Gclose(group);
  }

  return retval;
}

template <>
bool
H5ReadGroupAttribute<bool>(hid_t id,
                           const std::string& group_id,
                           const std::string& name,
                           bool& value)
{
  bool retval = false;

  auto group = H5Gopen2(id, group_id.c_str(), H5P_DEFAULT);
  if (group != H5I_INVALID_HID)
  {
    retval = H5ReadAttribute<bool>(group, name, value);
    H5Gclose(group);
  }

  return retval;
}

} // namespace opensn
