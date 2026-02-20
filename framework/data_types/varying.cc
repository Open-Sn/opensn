// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/data_types/varying.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <sstream>

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

void
Varying::CheckTypeMatch(const VaryingDataType type_A, const VaryingDataType type_B_required) const
{
  if (type_A != type_B_required)
    throw std::logic_error("Varying data type " + GetTypeName() +
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

Varying::Varying(const Varying& other) : type_(other.type_), data_(other.data_->Clone())
{
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
Varying::GetStringValue() const
{
  CheckTypeMatch(type_, VaryingDataType::STRING);

  return data_->GetStringValue();
}

bool
Varying::GetBoolValue() const
{
  CheckTypeMatch(type_, VaryingDataType::BOOL);

  return data_->GetBoolValue();
}

int64_t
Varying::GetIntegerValue() const
{
  CheckTypeMatch(type_, VaryingDataType::INTEGER);

  return data_->GetIntegerValue();
}

double
Varying::GetFloatValue() const
{
  CheckTypeMatch(type_, VaryingDataType::FLOAT);

  return data_->GetFloatValue();
}

size_t
Varying::GetByteSize() const
{
  return data_->Size();
}

std::string
Varying::PrintStr(bool with_type) const
{
  std::stringstream outstr;

  if (this->GetType() == VaryingDataType::STRING)
    outstr << "\"" << this->GetStringValue() << "\"";
  else if (this->GetType() == VaryingDataType::FLOAT)
    outstr << this->GetFloatValue() << (with_type ? "(double)" : "");
  else if (this->GetType() == VaryingDataType::INTEGER)
    outstr << this->GetIntegerValue();
  else if (this->GetType() == VaryingDataType::BOOL)
    outstr << (this->GetBoolValue() ? "true" : "false");

  return outstr.str();
}

} // namespace opensn

std::ostream&
operator<<(std::ostream& outstr, const opensn::Varying& value)
{
  outstr << value.PrintStr(false);
  return outstr;
}
