// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "H5Cpp.h"
#include <vector>
#include <optional>

namespace opensn
{
template <typename T>
H5::DataType get_datatype();

template <>
H5::DataType
get_datatype<char>()
{
  return H5::PredType::NATIVE_CHAR;
}

template <>
H5::DataType
get_datatype<signed char>()
{
  return H5::PredType::NATIVE_SCHAR;
}

template <>
H5::DataType
get_datatype<unsigned char>()
{
  return H5::PredType::NATIVE_UCHAR;
}

template <>
H5::DataType
get_datatype<short>()
{
  return H5::PredType::NATIVE_SHORT;
}

template <>
H5::DataType
get_datatype<unsigned short>()
{
  return H5::PredType::NATIVE_USHORT;
}

template <>
H5::DataType
get_datatype<int>()
{
  return H5::PredType::NATIVE_INT;
}

template <>
H5::DataType
get_datatype<unsigned int>()
{
  return H5::PredType::NATIVE_UINT;
}

template <>
H5::DataType
get_datatype<long>()
{
  return H5::PredType::NATIVE_LONG;
}

template <>
H5::DataType
get_datatype<unsigned long>()
{
  return H5::PredType::NATIVE_ULONG;
}

template <>
H5::DataType
get_datatype<long long>()
{
  return H5::PredType::NATIVE_LLONG;
}

template <>
H5::DataType
get_datatype<unsigned long long>()
{
  return H5::PredType::NATIVE_ULLONG;
}

template <>
H5::DataType
get_datatype<float>()
{
  return H5::PredType::NATIVE_FLOAT;
}

template <>
H5::DataType
get_datatype<double>()
{
  return H5::PredType::NATIVE_DOUBLE;
}

bool
H5Has(H5::H5File file, const std::string& name)
{
  return (H5Lexists(file.getId(), name.c_str(), H5P_DEFAULT) > 0);
}

template <typename T>
std::vector<T>
H5ReadDataset1D(H5::H5File file, const std::string& name)
{
  std::vector<T> data;

  try
  {
    H5::Exception::dontPrint();
    H5::DataSet dataset = file.openDataSet(name.c_str());
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);
    data.resize(dims[0]);
    dataset.read(data.data(), get_datatype<T>(), dataspace);
  }
  catch (H5::Exception error)
  {
    return std::vector<T>();
  }

  return data;
}

template <typename T>
bool
H5ReadAttribute(H5::H5File file, const std::string& name, T& value)
{
  try
  {
    H5::Attribute attribute = file.openAttribute(name.c_str());
    attribute.read(get_datatype<T>(), &value);
  }
  catch (H5::Exception error)
  {
    return false;
  }

  return true;
}

template <>
bool
H5ReadAttribute<std::string>(H5::H5File file, const std::string& name, std::string& value)
{
  try
  {
    H5::Attribute attribute = file.openAttribute(name.c_str());
    H5::StrType stype = attribute.getStrType();
    attribute.read(stype, value);
  }
  catch (H5::Exception error)
  {
    return false;
  }

  return true;
}

template <>
bool
H5ReadAttribute<bool>(H5::H5File file, const std::string& name, bool& value)
{
  try
  {
    H5::Attribute attribute = file.openAttribute(name.c_str());
    attribute.read(H5::PredType::NATIVE_INT, &value);
  }
  catch (H5::Exception error)
  {
    return false;
  }

  return true;
}

template <typename T>
bool
H5ReadGroupAttribute(H5::H5File file,
                     const std::string& group_id,
                     const std::string& name,
                     T& value)
{
  try
  {
    H5::Group group = file.openGroup(group_id.c_str());
    H5::Attribute attribute = group.openAttribute(name.c_str());
    attribute.read(get_datatype<T>(), &value);
  }
  catch (H5::Exception error)
  {
    return false;
  }

  return true;
}

template <>
bool
H5ReadGroupAttribute<std::string>(H5::H5File file,
                                  const std::string& group_id,
                                  const std::string& name,
                                  std::string& value)
{
  try
  {
    H5::Group group = file.openGroup(group_id.c_str());
    H5::Attribute attribute = group.openAttribute(name.c_str());
    H5::StrType stype = attribute.getStrType();
    attribute.read(stype, value);
  }
  catch (H5::Exception error)
  {
    return false;
  }

  return true;
}

template <>
bool
H5ReadGroupAttribute<bool>(H5::H5File file,
                           const std::string& group_id,
                           const std::string& name,
                           bool& value)
{
  try
  {
    H5::Group group = file.openGroup(group_id.c_str());
    H5::Attribute attribute = group.openAttribute(name.c_str());
    attribute.read(H5::PredType::NATIVE_INT, &value);
  }
  catch (H5::Exception error)
  {
    return false;
  }

  return true;
}
} // namespace opensn
