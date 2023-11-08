#include "framework/data_types/varying.h"

#include <algorithm>
#include <sstream>

#include "framework/logging/chi_log_exceptions.h"

/**Provides a string-name for an enumerated VaryingDataType.*/
std::string
chi_data_types::VaryingDataTypeStringName(chi_data_types::VaryingDataType type)
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

namespace chi_data_types
{
//=============================================== VaryingType Base
std::string
Varying::VaryingType::StringValue() const
{
  ChiLogicalError("Method not implemented");
}
bool
Varying::VaryingType::BoolValue() const
{
  ChiLogicalError("Method not implemented");
}
int64_t
Varying::VaryingType::IntegerValue() const
{
  ChiLogicalError("Method not implemented");
}
double
Varying::VaryingType::FloatValue() const
{
  ChiLogicalError("Method not implemented");
}
std::vector<std::byte>
Varying::VaryingType::BytesValue() const
{
  ChiLogicalError("Method not implemented");
}

//=============================================== VaryingByteArray
template <>
std::string
Varying::VaryingArbitraryType<std::vector<std::byte>>::StringValue() const
{
  ChiLogicalError("Method not implemented");
}
template <>
bool
Varying::VaryingArbitraryType<std::vector<std::byte>>::BoolValue() const
{
  ChiLogicalError("Method not implemented");
}
template <>
int64_t
Varying::VaryingArbitraryType<std::vector<std::byte>>::IntegerValue() const
{
  ChiLogicalError("Method not implemented");
}
template <>
double
Varying::VaryingArbitraryType<std::vector<std::byte>>::FloatValue() const
{
  ChiLogicalError("Method not implemented");
}

//=============================================== VaryingString
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
  ChiLogicalError("Method not implemented");
}
template <>
int64_t
Varying::VaryingArbitraryType<std::string>::IntegerValue() const
{
  ChiLogicalError("Method not implemented");
}
template <>
double
Varying::VaryingArbitraryType<std::string>::FloatValue() const
{
  ChiLogicalError("Method not implemented");
}

//=============================================== VaryingBool
template <>
std::string
Varying::VaryingArbitraryType<bool>::StringValue() const
{
  ChiLogicalError("Method not implemented");
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
  ChiLogicalError("Method not implemented");
}
template <>
double
Varying::VaryingArbitraryType<bool>::FloatValue() const
{
  ChiLogicalError("Method not implemented");
}

//=============================================== VaryingInteger
template <>
std::string
Varying::VaryingArbitraryType<int64_t>::StringValue() const
{
  ChiLogicalError("Method not implemented");
}
template <>
bool
Varying::VaryingArbitraryType<int64_t>::BoolValue() const
{
  ChiLogicalError("Method not implemented");
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
  ChiLogicalError("Method not implemented");
}

//=============================================== VaryingFloat
template <>
std::string
Varying::VaryingArbitraryType<double>::StringValue() const
{
  ChiLogicalError("Method not implemented");
}
template <>
bool
Varying::VaryingArbitraryType<double>::BoolValue() const
{
  ChiLogicalError("Method not implemented");
}
template <>
int64_t
Varying::VaryingArbitraryType<double>::IntegerValue() const
{
  ChiLogicalError("Method not implemented");
}
template <>
double
Varying::VaryingArbitraryType<double>::FloatValue() const
{
  return value_;
}

} // namespace chi_data_types

/**Checks if two VaryingDataType values match.
 * Type A is matched against type B.*/
void
chi_data_types::Varying::CheckTypeMatch(const VaryingDataType type_A,
                                        const VaryingDataType type_B_required) const
{
  if (type_A != type_B_required)
    throw std::logic_error("Varying data type " + TypeName() +
                           " does not "
                           "correspond to the required type, " +
                           VaryingDataTypeStringName(type_B_required));
}

// ###################################################################
// Constructors
/**Constructor for an arbitrary sequence of bytes value.*/
chi_data_types::Varying::Varying(const std::vector<std::byte>& value)
  : type_(VaryingDataType::ARBITRARY_BYTES)
{
  // raw_data_ = value;
  // data_initialized_ = true;
  data_ = std::make_unique<VaryingArbitraryType<std::vector<std::byte>>>(value);
}

/**Constructor for a string value.*/
chi_data_types::Varying::Varying(const std::string& value)
  : type_(VaryingDataType::STRING),
    data_(std::make_unique<VaryingArbitraryType<std::string>>(value))
{
}

/**Copy constructor.*/
chi_data_types::Varying::Varying(const Varying& other)
{
  data_ = other.data_->Clone();
  type_ = other.type_;
}

/**Move constructor.*/
chi_data_types::Varying::Varying(Varying&& other) noexcept
{
  std::swap(data_, other.data_);
  std::swap(type_, other.type_);
}

/**Assignment operator. i.e., type_A = type_B*/
chi_data_types::Varying&
chi_data_types::Varying::operator=(const Varying& other)
{
  if (this != &other)
  {
    data_ = other.data_->Clone();
    type_ = other.type_;
  }
  return *this;
}

// ###################################################################
//  Assignments
/**Assigns an arbitrary sequence of bytes value.*/
chi_data_types::Varying&
chi_data_types::Varying::operator=(const std::vector<std::byte>& value)
{
  type_ = VaryingDataType::ARBITRARY_BYTES;
  data_ = std::make_unique<VaryingArbitraryType<std::vector<std::byte>>>(value);
  return *this;
}

/**Assigns a string value.*/
chi_data_types::Varying&
chi_data_types::Varying::operator=(const std::string& value)
{
  type_ = VaryingDataType::STRING;
  data_ = std::make_unique<VaryingArbitraryType<std::string>>(value);
  return *this;
}

// ###################################################################
//  Get values
/**Returns the string value if valid. Otherwise throws std::logic_error.*/
std::string
chi_data_types::Varying::StringValue() const
{
  CheckTypeMatch(type_, VaryingDataType::STRING);

  return data_->StringValue();
}

/**Returns the bool value if valid. Otherwise throws std::logic_error.*/
bool
chi_data_types::Varying::BoolValue() const
{
  CheckTypeMatch(type_, VaryingDataType::BOOL);

  return data_->BoolValue();
}

/**Returns the integer value if valid. Otherwise throws std::logic_error.*/
int64_t
chi_data_types::Varying::IntegerValue() const
{
  CheckTypeMatch(type_, VaryingDataType::INTEGER);

  return data_->IntegerValue();
}

/**Returns the float value if valid. Otherwise throws std::logic_error.*/
double
chi_data_types::Varying::FloatValue() const
{
  CheckTypeMatch(type_, VaryingDataType::FLOAT);

  return data_->FloatValue();
}

// ###################################################################
/**Returns the raw byte size associated with the type.*/
size_t
chi_data_types::Varying::ByteSize() const
{
  return data_->Size();
}

// ###################################################################
/**Returns a string value for the value.*/
std::string
chi_data_types::Varying::PrintStr(bool with_type /*=false*/) const
{
  std::stringstream outstr;

  if (this->Type() == chi_data_types::VaryingDataType::STRING)
    outstr << "\"" << this->StringValue() << "\"";
  else if (this->Type() == chi_data_types::VaryingDataType::FLOAT)
    outstr << this->FloatValue() << (with_type ? "(double)" : "");
  else if (this->Type() == chi_data_types::VaryingDataType::INTEGER)
    outstr << this->IntegerValue();
  else if (this->Type() == chi_data_types::VaryingDataType::BOOL)
    outstr << (this->BoolValue() ? "true" : "false");

  return outstr.str();
}

// ###################################################################
/**Stream operator*/
std::ostream&
operator<<(std::ostream& outstr, const chi_data_types::Varying& value)
{
  outstr << value.PrintStr(false);
  return outstr;
}
