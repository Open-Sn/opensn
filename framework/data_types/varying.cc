// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/data_types/varying.h"

#include <algorithm>
#include <sstream>

#include "framework/logging/log_exceptions.h"

namespace opensn
{

std::string
VaryingDataTypeStringName(VaryingDataType type)
{
  switch (type)
  {
    case VaryingDataType::VOID:
      return "VOID";
    case VaryingDataType::ARBITRARY_BYTES:
      return "ARBITRARY_BYTES";
    case VaryingDataType::STRING:
      return "STRING";
    case VaryingDataType::BOOL:
      return "BOOL";
    case VaryingDataType::INTEGER:
      return "INTEGER";
    case VaryingDataType::FLOAT:
      return "FLOAT";
    default:
      return "UNKNOWN";
  }
}

// VaryingType Base
std::string
Varying::VaryingType::StringValue() const
{
  OpenSnLogicalError("Method not implemented");
}
bool
Varying::VaryingType::BoolValue() const
{
  OpenSnLogicalError("Method not implemented");
}
int64_t
Varying::VaryingType::IntegerValue() const
{
  OpenSnLogicalError("Method not implemented");
}
double
Varying::VaryingType::FloatValue() const
{
  OpenSnLogicalError("Method not implemented");
}
std::vector<std::byte>
Varying::VaryingType::BytesValue() const
{
  OpenSnLogicalError("Method not implemented");
}

// VaryingByteArray
template <>
std::string
Varying::VaryingArbitraryType<std::vector<std::byte>>::StringValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
bool
Varying::VaryingArbitraryType<std::vector<std::byte>>::BoolValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
int64_t
Varying::VaryingArbitraryType<std::vector<std::byte>>::IntegerValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
double
Varying::VaryingArbitraryType<std::vector<std::byte>>::FloatValue() const
{
  OpenSnLogicalError("Method not implemented");
}

// VaryingString
template <>
std::string
Varying::VaryingArbitraryType<std::string>::StringValue() const
{
  return value_;
}
template <>
bool
Varying::VaryingArbitraryType<std::string>::BoolValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
int64_t
Varying::VaryingArbitraryType<std::string>::IntegerValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
double
Varying::VaryingArbitraryType<std::string>::FloatValue() const
{
  OpenSnLogicalError("Method not implemented");
}

// VaryingBool
template <>
std::string
Varying::VaryingArbitraryType<bool>::StringValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
bool
Varying::VaryingArbitraryType<bool>::BoolValue() const
{
  return value_;
}
template <>
int64_t
Varying::VaryingArbitraryType<bool>::IntegerValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
double
Varying::VaryingArbitraryType<bool>::FloatValue() const
{
  OpenSnLogicalError("Method not implemented");
}

// VaryingInteger
template <>
std::string
Varying::VaryingArbitraryType<int64_t>::StringValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
bool
Varying::VaryingArbitraryType<int64_t>::BoolValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
int64_t
Varying::VaryingArbitraryType<int64_t>::IntegerValue() const
{
  return value_;
}
template <>
double
Varying::VaryingArbitraryType<int64_t>::FloatValue() const
{
  OpenSnLogicalError("Method not implemented");
}

// VaryingFloat
template <>
std::string
Varying::VaryingArbitraryType<double>::StringValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
bool
Varying::VaryingArbitraryType<double>::BoolValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
int64_t
Varying::VaryingArbitraryType<double>::IntegerValue() const
{
  OpenSnLogicalError("Method not implemented");
}
template <>
double
Varying::VaryingArbitraryType<double>::FloatValue() const
{
  return value_;
}

void
Varying::CheckTypeMatch(const VaryingDataType type_A, const VaryingDataType type_B_required) const
{
  if (type_A != type_B_required)
    throw std::logic_error("Varying data type " + TypeName() +
                           " does not "
                           "correspond to the required type, " +
                           VaryingDataTypeStringName(type_B_required));
}

// Constructors
Varying::Varying(const std::vector<std::byte>& value) : type_(VaryingDataType::ARBITRARY_BYTES)
{
  // raw_data_ = value;
  // data_initialized_ = true;
  data_ = std::make_unique<VaryingArbitraryType<std::vector<std::byte>>>(value);
}

Varying::Varying(const std::string& value)
  : type_(VaryingDataType::STRING),
    data_(std::make_unique<VaryingArbitraryType<std::string>>(value))
{
}

Varying::Varying(const Varying& other)
{
  data_ = other.data_->Clone();
  type_ = other.type_;
}

Varying::Varying(Varying&& other) noexcept
{
  std::swap(data_, other.data_);
  std::swap(type_, other.type_);
}

Varying&
Varying::operator=(const Varying& other)
{
  if (this != &other)
  {
    data_ = other.data_->Clone();
    type_ = other.type_;
  }
  return *this;
}

//  Assignments
Varying&
Varying::operator=(const std::vector<std::byte>& value)
{
  type_ = VaryingDataType::ARBITRARY_BYTES;
  data_ = std::make_unique<VaryingArbitraryType<std::vector<std::byte>>>(value);
  return *this;
}

Varying&
Varying::operator=(const std::string& value)
{
  type_ = VaryingDataType::STRING;
  data_ = std::make_unique<VaryingArbitraryType<std::string>>(value);
  return *this;
}

//  Get values
std::string
Varying::StringValue() const
{
  CheckTypeMatch(type_, VaryingDataType::STRING);

  return data_->StringValue();
}

bool
Varying::BoolValue() const
{
  CheckTypeMatch(type_, VaryingDataType::BOOL);

  return data_->BoolValue();
}

int64_t
Varying::IntegerValue() const
{
  CheckTypeMatch(type_, VaryingDataType::INTEGER);

  return data_->IntegerValue();
}

double
Varying::FloatValue() const
{
  CheckTypeMatch(type_, VaryingDataType::FLOAT);

  return data_->FloatValue();
}

size_t
Varying::ByteSize() const
{
  return data_->Size();
}

std::string
Varying::PrintStr(bool with_type) const
{
  std::stringstream outstr;

  if (this->Type() == VaryingDataType::STRING)
    outstr << "\"" << this->StringValue() << "\"";
  else if (this->Type() == VaryingDataType::FLOAT)
    outstr << this->FloatValue() << (with_type ? "(double)" : "");
  else if (this->Type() == VaryingDataType::INTEGER)
    outstr << this->IntegerValue();
  else if (this->Type() == VaryingDataType::BOOL)
    outstr << (this->BoolValue() ? "true" : "false");

  return outstr.str();
}

} // namespace opensn

std::ostream&
operator<<(std::ostream& outstr, const opensn::Varying& value)
{
  outstr << value.PrintStr(false);
  return outstr;
}
