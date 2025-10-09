// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include "framework/logging/log_exceptions.h"
#include <cstddef>
#include <iostream>
#include <type_traits>
#include <vector>
#include <memory>

namespace opensn
{

/// Enumeration of data-types supported by Varying
enum class VaryingDataType : int
{
  VOID = 0,            ///< Basically undefined or null
  ARBITRARY_BYTES = 1, ///< Basic sequence of bytes
  STRING = 2,          ///< Datatype mapping to std::string
  BOOL = 3,            ///< Datatype mapping to bool
  INTEGER = 4,         ///< Datatype mapping to int64_t
  FLOAT = 5,           ///< Datatype mapping to double
  USER_DATA = 6        ///< Datatype mapping to user data
};

/// Provides a string-name for an enumerated VaryingDataType.
std::string VaryingDataTypeStringName(VaryingDataType type);

template <typename T>
struct is_shared_ptr : std::false_type
{
};

template <typename T>
struct is_shared_ptr<std::shared_ptr<T>> : std::true_type
{
};

template <typename T>
inline constexpr bool is_shared_ptr_v = is_shared_ptr<T>::value;

class Varying
{
private:
  // Helpers
  template <typename T>
  struct IsByteArray
  {
    static constexpr bool value = std::is_same_v<T, std::vector<std::byte>>;
  };
  template <typename T>
  struct IsBool
  {
    static constexpr bool value = std::is_same_v<T, bool>;
  };
  template <typename T>
  struct IsString
  {
    static constexpr bool value = std::is_same_v<T, std::string> or std::is_same_v<T, char*>;
  };
  template <typename T>
  struct IsFloat
  {
    static constexpr bool value = std::is_floating_point_v<T>;
  };
  template <typename T>
  struct IsInteger
  {
    static constexpr bool value = std::is_integral_v<T> and not std::is_same_v<T, bool>;
  };
  template <typename T>
  struct IsUserData
  {
    static constexpr bool value = std::is_pointer_v<T> or is_shared_ptr_v<T> or
                                  (std::is_class_v<T> and not std::is_same_v<T, std::string>);
  };

  template <typename T>
  using BoolType = typename std::enable_if_t<IsBool<T>::value, T>;
  template <typename T>
  using FloatType = typename std::enable_if_t<IsFloat<T>::value, T>;
  template <typename T>
  using IntegerType = typename std::enable_if_t<IsInteger<T>::value, T>;

  template <typename T>
  using BoolStorageType = typename std::enable_if_t<IsBool<T>::value, bool>;
  template <typename T>
  using FloatStorageType = typename std::enable_if_t<IsFloat<T>::value, double>;
  template <typename T>
  using IntegerStorageType = typename std::enable_if_t<IsInteger<T>::value, int64_t>;
  template <typename T>
  using UserDataStorageType = typename std::enable_if_t<IsUserData<T>::value, T>;

  template <typename T>
  BoolStorageType<T> CastValue(const T& value)
  {
    return value;
  }

  template <typename T>
  FloatStorageType<T> CastValue(const T& value)
  {
    return static_cast<double>(value);
  }

  template <typename T>
  IntegerStorageType<T> CastValue(const T& value)
  {
    return static_cast<int64_t>(value);
  }

  template <typename T>
  UserDataStorageType<T> CastValue(const T& value)
  {
    return static_cast<T>(value);
  }

  /// This acts as a base class for templated child arbitrary types
  class VaryingType
  {
  public:
    virtual std::string GetStringValue() const = 0;
    virtual bool GetBoolValue() const = 0;
    virtual int64_t GetIntegerValue() const = 0;
    virtual double GetFloatValue() const = 0;
    virtual std::vector<std::byte> GetBytesValue() const = 0;

    virtual std::unique_ptr<VaryingType> Clone() const = 0;
    virtual size_t Size() const = 0;

    virtual bool operator==(const VaryingType& that) const = 0;
    virtual bool operator!=(const VaryingType& that) const = 0;
    virtual bool operator>(const VaryingType& that) const = 0;
    virtual bool operator<(const VaryingType& that) const = 0;
    virtual bool operator>=(const VaryingType& that) const = 0;
    virtual bool operator<=(const VaryingType& that) const = 0;

    VaryingDataType GetType() const { return type_; }

    virtual ~VaryingType() = default;

  protected:
    VaryingDataType type_;
    explicit VaryingType(VaryingDataType type) : type_(type) {}
  };

  template <typename T>
  class VaryingArbitraryType : public VaryingType
  {
  public:
    explicit VaryingArbitraryType(T value)
      : VaryingType(IsByteArray<T>::value  ? VaryingDataType::ARBITRARY_BYTES
                    : IsString<T>::value   ? VaryingDataType::STRING
                    : IsBool<T>::value     ? VaryingDataType::BOOL
                    : IsInteger<T>::value  ? VaryingDataType::INTEGER
                    : IsFloat<T>::value    ? VaryingDataType::FLOAT
                    : IsUserData<T>::value ? VaryingDataType::USER_DATA
                                           : VaryingDataType::VOID),
        value_(std::move(value))
    {
    }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    std::string GetStringValue() const override { OpenSnLogicalError("Method not implemented"); }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    bool GetBoolValue() const override { OpenSnLogicalError("Method not implemented"); }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    int64_t GetIntegerValue() const override { OpenSnLogicalError("Method not implemented"); }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    double GetFloatValue() const override { OpenSnLogicalError("Method not implemented"); }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    std::vector<std::byte> GetBytesValue() const override
    {
      OpenSnLogicalError("Method not implemented");
    }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    std::unique_ptr<VaryingType> Clone() const override
    {
      return std::make_unique<VaryingArbitraryType<T>>(value_);
    }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    size_t Size() const override { return sizeof(T); }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    bool operator==(const VaryingType& that) const override
    {
      if (type_ != that.GetType())
        return false;

      switch (this->GetType())
      {
        case VaryingDataType::ARBITRARY_BYTES:
          return GetBytesValue() == that.GetBytesValue();
        case VaryingDataType::STRING:
          return GetStringValue() == that.GetStringValue();
        case VaryingDataType::BOOL:
          return GetBoolValue() == that.GetBoolValue();
        case VaryingDataType::INTEGER:
          return GetIntegerValue() == that.GetIntegerValue();
        case VaryingDataType::FLOAT:
          return GetFloatValue() == that.GetFloatValue();
        case VaryingDataType::VOID:
        default:
          return false;
      }
    }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    bool operator!=(const VaryingType& that) const override { return not(*this == that); }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    bool operator>(const VaryingType& that) const override
    {
      if (type_ != that.GetType())
        return false;

      switch (this->GetType())
      {
        case VaryingDataType::ARBITRARY_BYTES:
          return GetBytesValue() > that.GetBytesValue();
        case VaryingDataType::STRING:
          return GetStringValue() > that.GetStringValue();
        case VaryingDataType::BOOL:
          return GetBoolValue() > that.GetBoolValue();
        case VaryingDataType::INTEGER:
          return GetIntegerValue() > that.GetIntegerValue();
        case VaryingDataType::FLOAT:
          return GetFloatValue() > that.GetFloatValue();
        case VaryingDataType::VOID:
        default:
          return false;
      }
    }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    bool operator<(const VaryingType& that) const override
    {
      if (type_ != that.GetType())
        return false;

      switch (this->GetType())
      {
        case VaryingDataType::ARBITRARY_BYTES:
          return GetBytesValue() < that.GetBytesValue();
        case VaryingDataType::STRING:
          return GetStringValue() < that.GetStringValue();
        case VaryingDataType::BOOL:
          return GetBoolValue() < that.GetBoolValue();
        case VaryingDataType::INTEGER:
          return GetIntegerValue() < that.GetIntegerValue();
        case VaryingDataType::FLOAT:
          return GetFloatValue() < that.GetFloatValue();
        case VaryingDataType::VOID:
        default:
          return false;
      }
    }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    bool operator>=(const VaryingType& that) const override
    {
      return (*this > that) or (*this == that);
    }

    // NOLINTNEXTLINE(portability-template-virtual-member-function)
    bool operator<=(const VaryingType& that) const override
    {
      return (*this < that) or (*this == that);
    }

    T Value() const { return value_; }

  private:
    T value_;
  };

  /// Type specification
  VaryingDataType type_ = VaryingDataType::VOID;
  std::unique_ptr<VaryingType> data_ = nullptr;

private:
  /**
   * Checks if two VaryingDataType values match.
   * Type A is matched against type B.
   */
  void CheckTypeMatch(VaryingDataType type_A, VaryingDataType type_B_required) const;

private:
public:
  // Constructors
  /**
   * Generalized constructor for bool, integral- and float-types. This constructor has been
   * specialized for std::string and std::vector<std::byte>.
   */
  template <typename T>
  explicit Varying(const T& value)
  {
    constexpr bool is_supported_type =
      IsBool<T>::value or IsFloat<T>::value or IsInteger<T>::value or IsUserData<T>::value;
    static_assert(is_supported_type, "Constructor called with unsupported type");

    if (IsBool<T>::value)
    {
      type_ = VaryingDataType::BOOL;
    }
    else if (IsFloat<T>::value)
    {
      type_ = VaryingDataType::FLOAT;
    }
    else if (IsInteger<T>::value)
    {
      type_ = VaryingDataType::INTEGER;
    }
    else
      type_ = VaryingDataType::USER_DATA;

    data_ = Helper(CastValue(value));
  }

  static std::unique_ptr<VaryingType> Helper(const bool& value)
  {
    return std::make_unique<VaryingArbitraryType<bool>>(value);
  }

  static std::unique_ptr<VaryingType> Helper(const int64_t& value)
  {
    return std::make_unique<VaryingArbitraryType<int64_t>>(value);
  }

  static std::unique_ptr<VaryingType> Helper(const double& value)
  {
    return std::make_unique<VaryingArbitraryType<double>>(value);
  }

  static std::unique_ptr<VaryingType> Helper(const Vector3& value)
  {
    return std::make_unique<VaryingArbitraryType<Vector3>>(value);
  }

  template <typename T>
  static std::unique_ptr<VaryingType> Helper(const std::shared_ptr<T>& value)
  {
    return std::make_unique<VaryingArbitraryType<std::shared_ptr<T>>>(value);
  }

  /// Constructor for an arbitrary sequence of bytes value.
  explicit Varying(const std::vector<std::byte>& value);

  /// Constructor for a string value.
  explicit Varying(const std::string& value);

  /// Constructor for a string literal value.
  explicit Varying(const char* value) : Varying((not value) ? std::string() : std::string(value)) {}

  /// Copy constructor.
  Varying(const Varying& other);

  /// Move constructor.
  Varying(Varying&& other) noexcept;

public:
  // Copy assignment operator
  /// Assignment operator. i.e., type_A = type_B
  Varying& operator=(const Varying& other);

  // Assignment operators
  /// Assigns an arbitrary sequence of bytes value.
  Varying& operator=(const std::vector<std::byte>& value);

  /// Assigns a string value.
  Varying& operator=(const std::string& value);

  /// Assigns a bool value.
  template <typename T, std::enable_if_t<IsBool<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::BOOL;
    data_ = std::make_unique<VaryingArbitraryType<bool>>(value);

    return *this;
  }

  /// Assigns an integer value.
  template <typename T, std::enable_if_t<IsInteger<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::INTEGER;
    data_ = std::make_unique<VaryingArbitraryType<int64_t>>(value);

    return *this;
  }

  /// Assign a floating point value.
  template <typename T, std::enable_if_t<IsFloat<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::FLOAT;
    data_ = std::make_unique<VaryingArbitraryType<double>>(value);

    return *this;
  }

  /// Equality operator
  bool operator==(const Varying& that) const { return *data_ == *that.data_; }

  /// Inequality operator
  bool operator!=(const Varying& that) const { return not(*this == that); }

  /// Relation operators
  bool operator>(const Varying& that) const { return *data_ > *that.data_; }

  /// Relation operators
  bool operator>=(const Varying& that) const { return (*this > that) or (*this == that); }

  /// Relation operators
  bool operator<(const Varying& that) const { return *data_ < *that.data_; }

  /// Relation operators
  bool operator<=(const Varying& that) const { return (*this < that) or (*this == that); }

  /// Returns a default value for the type required.
  template <typename T>
  static T DefaultValue()
  {
    return {};
  }

public:
  // More Helpers
  template <typename T>
  struct IsSignedInteger
  {
    static constexpr bool value =
      std::is_integral_v<T> and std::is_signed_v<T> and not std::is_same_v<T, bool>;
  };
  template <typename T>
  struct IsUnsignedInteger
  {
    static constexpr bool value =
      std::is_integral_v<T> and std::is_unsigned_v<T> and not std::is_same_v<T, bool>;
  };

  template <typename T>
  using StringType = typename std::enable_if_t<IsString<T>::value, T>;
  template <typename T>
  using SignedIntegerType = typename std::enable_if_t<IsSignedInteger<T>::value, T>;
  template <typename T>
  using UnsignedIntegerType = typename std::enable_if_t<IsUnsignedInteger<T>::value, T>;
  template <typename T>
  using UserDataType = typename std::enable_if_t<IsUserData<T>::value, T>;

  /// Returns values of type bool if able.
  template <typename T>
  BoolType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::BOOL);

    return data_->GetBoolValue();
  }

  /// Returns floating point values if able.
  template <typename T>
  FloatType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::FLOAT);

    const double value = data_->GetFloatValue();

    return static_cast<T>(value);
  }

  /// Returns a string if able.
  template <typename T>
  StringType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::STRING);

    return data_->GetStringValue();
  }

  /// Returns a signed integer if able.
  template <typename T>
  SignedIntegerType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::INTEGER);

    const int64_t value = data_->GetIntegerValue();

    return static_cast<T>(value);
  }

  /// Returns an unsigned integer if able.
  template <typename T>
  UnsignedIntegerType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::INTEGER);

    const int64_t value = data_->GetIntegerValue();

    if (value < 0)
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Attempt to cast negative number to unsigned.");

    return static_cast<T>(value);
  }

  /// Returns user data if able.
  template <typename T>
  UserDataType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::USER_DATA);

    const auto value = dynamic_cast<VaryingArbitraryType<T>*>(data_.get())->Value();

    return static_cast<T>(value);
  }

  /// Returns the string value if valid. Otherwise throws std::logic_error.
  std::string GetStringValue() const;

  /// Returns the bool value if valid. Otherwise throws std::logic_error.
  bool GetBoolValue() const;

  /// Returns the integer value if valid. Otherwise throws std::logic_error.
  int64_t GetIntegerValue() const;

  /// Returns the float value if valid. Otherwise throws std::logic_error.
  double GetFloatValue() const;

  template <typename T>
  T GetUserDataValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::USER_DATA);
    return GetValue<T>();
  }

  /// Returns the raw byte size associated with the type.
  size_t GetByteSize() const;

public:
  /// Returns the current-type of the variable.
  VaryingDataType GetType() const { return type_; }

  /// Returns the string type name of the type.
  std::string GetTypeName() const { return VaryingDataTypeStringName(type_); }

  /// Returns a string value for the value.
  std::string PrintStr(bool with_type = true) const;

public:
  ~Varying() = default;
};

template <>
inline std::string
Varying::VaryingArbitraryType<std::string>::GetStringValue() const
{
  return value_;
}

template <>
inline bool
Varying::VaryingArbitraryType<bool>::GetBoolValue() const
{
  return value_;
}

template <>
inline int64_t
Varying::VaryingArbitraryType<int64_t>::GetIntegerValue() const
{
  return value_;
}

template <>
inline double
Varying::VaryingArbitraryType<double>::GetFloatValue() const
{
  return value_;
}

} // namespace opensn

/// Stream operator
std::ostream& operator<<(std::ostream& outstr, const opensn::Varying& value);
