#pragma once

#include "opensn/framework/data_types/chi_data_types.h"

#include <cstddef>
#include <iostream>
#include <vector>
#include <memory>

namespace chi_data_types
{
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

  /**This acts as a base class for templated child arbitrary types*/
  class VaryingType
  {
  public:
    virtual std::string StringValue() const;
    virtual bool BoolValue() const;
    virtual int64_t IntegerValue() const;
    virtual double FloatValue() const;
    virtual std::vector<std::byte> BytesValue() const;

    virtual std::unique_ptr<VaryingType> Clone() const = 0;
    virtual size_t Size() const = 0;

    virtual bool operator==(const VaryingType& that) const = 0;
    virtual bool operator!=(const VaryingType& that) const = 0;
    virtual bool operator>(const VaryingType& that) const = 0;
    virtual bool operator<(const VaryingType& that) const = 0;
    virtual bool operator>=(const VaryingType& that) const = 0;
    virtual bool operator<=(const VaryingType& that) const = 0;

    VaryingDataType Type() const { return type_; }

    virtual ~VaryingType() = default;

  protected:
    VaryingDataType type_;
    explicit VaryingType(VaryingDataType type) : type_(type) {}
  };

  template <typename T>
  class VaryingArbitraryType : public VaryingType
  {
  public:
    // clang-format off
    explicit VaryingArbitraryType(T value)
      : VaryingType(IsByteArray<T>::value ? VaryingDataType::ARBITRARY_BYTES :
                    IsString<T>::value ? VaryingDataType::STRING :
                    IsBool<T>::value ? VaryingDataType::BOOL :
                    IsInteger<T>::value ? VaryingDataType::INTEGER :
                    IsFloat<T>::value ? VaryingDataType::FLOAT :
                    VaryingDataType::VOID),
      value_(value)
    {
    }
    // clang-format on
    std::string StringValue() const override;
    bool BoolValue() const override;
    int64_t IntegerValue() const override;
    double FloatValue() const override;

    std::unique_ptr<VaryingType> Clone() const override
    {
      return std::make_unique<VaryingArbitraryType<T>>(value_);
    }
    size_t Size() const override { return sizeof(T); }

    bool operator==(const VaryingType& that) const override
    {
      if (type_ != that.Type()) return false;

      switch (this->Type())
      {
        case VaryingDataType::ARBITRARY_BYTES:
          return BytesValue() == that.BytesValue();
        case VaryingDataType::STRING:
          return StringValue() == that.StringValue();
        case VaryingDataType::BOOL:
          return BoolValue() == that.BoolValue();
        case VaryingDataType::INTEGER:
          return IntegerValue() == that.IntegerValue();
        case VaryingDataType::FLOAT:
          return FloatValue() == that.FloatValue();
        case VaryingDataType::VOID:
        default:
          return false;
      }
    }

    bool operator!=(const VaryingType& that) const override { return not(*this == that); }
    bool operator>(const VaryingType& that) const override
    {
      if (type_ != that.Type()) return false;

      switch (this->Type())
      {
        case VaryingDataType::ARBITRARY_BYTES:
          return BytesValue() > that.BytesValue();
        case VaryingDataType::STRING:
          return StringValue() > that.StringValue();
        case VaryingDataType::BOOL:
          return BoolValue() > that.BoolValue();
        case VaryingDataType::INTEGER:
          return IntegerValue() > that.IntegerValue();
        case VaryingDataType::FLOAT:
          return FloatValue() > that.FloatValue();
        case VaryingDataType::VOID:
        default:
          return false;
      }
    }
    bool operator<(const VaryingType& that) const override
    {
      if (type_ != that.Type()) return false;

      switch (this->Type())
      {
        case VaryingDataType::ARBITRARY_BYTES:
          return BytesValue() < that.BytesValue();
        case VaryingDataType::STRING:
          return StringValue() < that.StringValue();
        case VaryingDataType::BOOL:
          return BoolValue() < that.BoolValue();
        case VaryingDataType::INTEGER:
          return IntegerValue() < that.IntegerValue();
        case VaryingDataType::FLOAT:
          return FloatValue() < that.FloatValue();
        case VaryingDataType::VOID:
        default:
          return false;
      }
    }
    bool operator>=(const VaryingType& that) const override
    {
      return (*this > that) or (*this == that);
    }
    bool operator<=(const VaryingType& that) const override
    {
      return (*this < that) or (*this == that);
    }

  private:
    T value_;
  };

  /**Type specification*/
  VaryingDataType type_ = VaryingDataType::VOID;
  std::unique_ptr<VaryingType> data_ = nullptr;

private:
  /**Checks if two VaryingDataType values match.
   * Type A is matched against type B.*/
  void CheckTypeMatch(VaryingDataType type_A, VaryingDataType type_B_required) const;

private:
public:
  // Constructors
  /**Generalized constructor for bool, integral- and float-types. This
   * constructor has been specialized for std::string and
   * std::vector<std::byte>.*/
  template <typename T>
  explicit Varying(const T& value)
  {
    constexpr bool is_supported_type = IsBool<T>::value or IsFloat<T>::value or IsInteger<T>::value;
    static_assert(is_supported_type, "Constructor called with unsupported type");

    if (IsBool<T>::value) { type_ = VaryingDataType::BOOL; }
    else if (IsFloat<T>::value) { type_ = VaryingDataType::FLOAT; }
    else if (IsInteger<T>::value) { type_ = VaryingDataType::INTEGER; }

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

  // template <typename T>
  // explicit Varying(const BoolType<T>& value)
  //{
  //   constexpr bool is_supported_type = IsBool<T>::value;
  //   static_assert(is_supported_type,
  //                 "Constructor called with unsupported type");
  //
  //   type_ = VaryingDataType::BOOL;
  //   data_ = std::make_unique<VaryingArbitraryType<BoolStorageType<T>>>(
  //     CastValue(value));
  // }
  //
  // template <typename T>
  // explicit Varying(const IntegerType<T>& value)
  //{
  //   constexpr bool is_supported_type = IsInteger<T>::value;
  //   static_assert(is_supported_type,
  //                 "Constructor called with unsupported type");
  //
  //   type_ = VaryingDataType::INTEGER;
  //   data_ = std::make_unique<VaryingArbitraryType<IntegerStorageType<T>>>(
  //     CastValue(value));
  // }
  //
  // template <typename T>
  // explicit Varying(const FloatType<T>& value)
  //{
  //   constexpr bool is_supported_type = IsFloat<T>::value;
  //   static_assert(is_supported_type,
  //                 "Constructor called with unsupported type");
  //
  //   type_ = VaryingDataType::FLOAT;
  //   data_ = std::make_unique<VaryingArbitraryType<FloatStorageType<T>>>(
  //     CastValue(value));
  // }

  /**Constructor for an arbitrary sequence of bytes value.*/
  explicit Varying(const std::vector<std::byte>& value);
  /**Constructor for a string value.*/
  explicit Varying(const std::string& value);
  /**Constructor for a string literal value.*/
  // explicit Varying(const char*& value) : Varying(std::string(value)) {}
  explicit Varying(const char* value) : Varying((not value) ? std::string() : std::string(value)) {}
  // template <std::size_t N>
  // explicit Varying(const char (&value)[N])
  //   : Varying(static_cast<const char*>(value))
  //{
  // }

  /**Copy constructor.*/
  Varying(const Varying& other);

  /**Move constructor.*/
  Varying(Varying&& other) noexcept;

public:
  // Copy assignment operator
  /**Assignment operator. i.e., type_A = type_B*/
  Varying& operator=(const Varying& other);

  // Assignment operators
  /**Assigns an arbitrary sequence of bytes value.*/
  Varying& operator=(const std::vector<std::byte>& value);
  /**Assigns a string value.*/
  Varying& operator=(const std::string& value);

  /**Assigns a bool value.*/
  template <typename T, std::enable_if_t<IsBool<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::BOOL;
    data_ = std::make_unique<VaryingArbitraryType<bool>>(value);

    return *this;
  }

  /**Assigns an integer value.*/
  template <typename T, std::enable_if_t<IsInteger<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::INTEGER;
    data_ = std::make_unique<VaryingArbitraryType<int64_t>>(value);

    return *this;
  }

  /**Assign a floating point value.*/
  template <typename T, std::enable_if_t<IsFloat<T>::value, bool> = true>
  Varying& operator=(const T& value)
  {
    type_ = VaryingDataType::FLOAT;
    data_ = std::make_unique<VaryingArbitraryType<double>>(value);

    return *this;
  }

  /**Equality operator*/
  bool operator==(const Varying& that) const { return *data_ == *that.data_; }

  /**Inequality operator*/
  bool operator!=(const Varying& that) const { return not(*this == that); }

  /**Relation operators*/
  bool operator>(const Varying& that) const { return *data_ > *that.data_; }
  /**Relation operators*/
  bool operator>=(const Varying& that) const { return (*this > that) or (*this == that); }
  /**Relation operators*/
  bool operator<(const Varying& that) const { return *data_ < *that.data_; }
  /**Relation operators*/
  bool operator<=(const Varying& that) const { return (*this < that) or (*this == that); }

  /**Returns a default value for the type required.*/
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

  /**Returns values of type bool if able.*/
  template <typename T>
  BoolType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::BOOL);

    return data_->BoolValue();
  }

  /**Returns floating point values if able.*/
  template <typename T>
  FloatType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::FLOAT);

    const double value = data_->FloatValue();

    return static_cast<T>(value);
  }

  /**Returns a string if able.*/
  template <typename T>
  StringType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::STRING);

    return data_->StringValue();
  }

  /**Returns a signed integer if able.*/
  template <typename T>
  SignedIntegerType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::INTEGER);

    const int64_t value = data_->IntegerValue();

    return static_cast<T>(value);
  }

  /**Returns an unsigned integer if able.*/
  template <typename T>
  UnsignedIntegerType<T> GetValue() const
  {
    CheckTypeMatch(type_, VaryingDataType::INTEGER);

    const int64_t value = data_->IntegerValue();

    if (value < 0)
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Attempt to cast negative number to unsigned.");

    return static_cast<T>(value);
  }

  /**Returns the string value if valid. Otherwise throws std::logic_error.*/
  std::string StringValue() const;
  /**Returns the bool value if valid. Otherwise throws std::logic_error.*/
  bool BoolValue() const;
  /**Returns the integer value if valid. Otherwise throws std::logic_error.*/
  int64_t IntegerValue() const;
  /**Returns the float value if valid. Otherwise throws std::logic_error.*/
  double FloatValue() const;

  /**Returns the raw byte size associated with the type.*/
  size_t ByteSize() const;

public:
  /**Returns the current-type of the variable.*/
  VaryingDataType Type() const { return type_; }
  /**Returns the string type name of the type.*/
  std::string TypeName() const { return VaryingDataTypeStringName(type_); }

  /**Returns a string value for the value.*/
  std::string PrintStr(bool with_type = true) const;

public:
  ~Varying() = default;
}; // class Varying

} // namespace chi_data_types

/**Stream operator*/
std::ostream& operator<<(std::ostream& outstr, const chi_data_types::Varying& value);
