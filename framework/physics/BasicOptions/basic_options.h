#ifndef CHI_BASIC_OPTIONS_H
#define CHI_BASIC_OPTIONS_H

#include "data_types/varying.h"

namespace chi_physics
{
//###################################################################
/**Class for option.*/
class BasicOption
{
private:
  std::string             name_;
  chi_data_types::Varying value_;

public:
  BasicOption(const std::string& name, const std::string& string_value) :
    name_(name), value_(string_value) {}

  BasicOption(const std::string& name, const bool& bool_value) :
    name_(name), value_(bool_value) {}

  BasicOption(const std::string& name, const int64_t& integer_value) :
    name_(name), value_(integer_value) {}

  BasicOption(const std::string& name, const double& float_value) :
    name_(name), value_(float_value) {}

  chi_data_types::VaryingDataType Type() const {return value_.Type();}

  std::string Name() const {return name_;}
  std::string StringValue() const  {return value_.StringValue();}
  bool        BoolValue() const    {return value_.BoolValue();}
  int64_t     IntegerValue() const {return value_.IntegerValue();}
  double      FloatValue() const   {return value_.FloatValue();}

  void SetStringValue(const std::string& in_string_value) { value_ = in_string_value;}
  void SetBoolValue(const bool& in_bool_value)            { value_ = in_bool_value;}
  void SetIntegerValue(const int64_t& in_integer_value)   { value_ = in_integer_value;}
  void SetFloatValue(const double& in_float_value)        { value_ = in_float_value;}

};

//###################################################################
/**Class for basic options*/
class BasicOptions
{
private:
  std::vector<BasicOption> options_;

public:
  BasicOptions() = default;

  /**Constructor with initializer list.*/
  BasicOptions(std::initializer_list<BasicOption> in_options) :
    options_(in_options)
  {  }

  //Operators
  /**Returns a constant reference to an option that matches the
   * requested name. If no name-match is found the method will throw
   * a std::out_of_range exception.*/
  const BasicOption& operator()(const std::string& option_name) const;

  /**Returns a constant reference to an option at the given
   * index. If the index is out of range then a std::out_of_range
   * exception is thrown. This method can potentially be faster than
   * the string comparison equivalent.*/
  const BasicOption& operator()(size_t index) const;

  /**Returns a non-constant reference to an option that matches the
   * requested name. If no name-match is found the method will throw
   * a std::out_of_range exception.*/
  BasicOption& operator[](const std::string& option_name);

  /**Returns a non-constant reference to an option at the given
   * index. If the index is out of range then a std::out_of_range
   * exception is thrown. This method can potentially be faster than
   * the string comparison equivalent.*/
  BasicOption& operator[](size_t index);

  //AddOption
  template<typename T>
  void AddOption(const std::string& option_name, const T& value);

  //Utilities
  /**Attempts to find an option that matches the requested name.
   * If one is found then its corresponding index is
   * returned. If it is not found then a std::out_of_range
   * exception is thrown.*/
  size_t GetOptionIndexFromName(const std::string& option_name) const;

};

}//namespace chi_physics
#endif