// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <utility>

#include "framework/data_types/varying.h"

namespace opensn
{

/// Class for option.
class BasicOption
{
private:
  std::string name_;
  Varying value_;

public:
  BasicOption(std::string name, const std::string& string_value)
    : name_(std::move(name)), value_(string_value)
  {
  }

  BasicOption(std::string name, const bool& bool_value) : name_(std::move(name)), value_(bool_value)
  {
  }

  BasicOption(std::string name, const int64_t& integer_value)
    : name_(std::move(name)), value_(integer_value)
  {
  }

  BasicOption(std::string name, const double& float_value)
    : name_(std::move(name)), value_(float_value)
  {
  }

  VaryingDataType Type() const { return value_.Type(); }

  std::string Name() const { return name_; }
  std::string StringValue() const { return value_.StringValue(); }
  bool BoolValue() const { return value_.BoolValue(); }
  int64_t IntegerValue() const { return value_.IntegerValue(); }
  double FloatValue() const { return value_.FloatValue(); }

  void SetStringValue(const std::string& value) { value_ = value; }
  void SetBoolValue(const bool& value) { value_ = value; }
  void SetIntegerValue(const int64_t& value) { value_ = value; }
  void SetFloatValue(const double& value) { value_ = value; }
};

/// Class for basic options
class BasicOptions
{
private:
  std::vector<BasicOption> options_;

public:
  BasicOptions() = default;

  /// Constructor with initializer list.
  BasicOptions(std::initializer_list<BasicOption> options) : options_(options) {}

  // Operators

  /**
   * Returns a constant reference to an option that matches the requested name. If no name-match is
   * found the method will throw a std::out_of_range exception.
   */
  const BasicOption& operator()(const std::string& option_name) const;

  /**
   * Returns a constant reference to an option at the given index. If the index is out of range then
   * a std::out_of_range exception is thrown. This method can potentially be faster than the string
   * comparison equivalent.
   */
  const BasicOption& operator()(size_t index) const;

  /**
   * Returns a non-constant reference to an option that matches the requested name. If no name-match
   * is found the method will throw a std::out_of_range exception.
   */
  BasicOption& operator[](const std::string& option_name);

  /**
   * Returns a non-constant reference to an option at the given index. If the index is out of range
   * then a std::out_of_range exception is thrown. This method can potentially be faster than the
   * string comparison equivalent.
   */
  BasicOption& operator[](size_t index);

  /// Adds an option to the options list.
  template <typename T>
  void AddOption(const std::string& option_name, const T& value);

  // Utilities

  /**
   * Attempts to find an option that matches the requested name. If one is found then its
   * corresponding index is returned. If it is not found then a std::out_of_range exception is
   * thrown.
   */
  size_t GetOptionIndexFromName(const std::string& option_name) const;
};

} // namespace opensn
