// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/parameters/parameter_block.h"
#include <algorithm>
#include <memory>
#include <sstream>
#include <iostream>

namespace opensn
{

std::string
ParameterBlockTypeName(ParameterBlockType type)
{
  switch (type)
  {
    case ParameterBlockType::BOOLEAN:
      return "BOOLEAN";
    case ParameterBlockType::FLOAT:
      return "FLOAT";
    case ParameterBlockType::STRING:
      return "STRING";
    case ParameterBlockType::INTEGER:
      return "INTEGER";
    case ParameterBlockType::ARRAY:
      return "ARRAY";
    case ParameterBlockType::BLOCK:
      return "BLOCK";
    case ParameterBlockType::USER_DATA:
      return "USER_DATA";
    default:
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) + ": No name associated with type");
  }
}

void
ParameterBlock::SetBlockName(const std::string& name)
{
  name_ = name;
}

ParameterBlock::ParameterBlock(const std::string& name)
  : type_(ParameterBlockType::BLOCK), name_(name), value_ptr_(nullptr)
{
}

ParameterBlock::ParameterBlock(const ParameterBlock& other)
  : type_(other.type_),
    name_(other.name_),
    parameters_(other.parameters_),
    error_origin_scope_(other.error_origin_scope_)
{
  if (other.value_ptr_)
    value_ptr_ = std::make_unique<Varying>(*other.value_ptr_);
}

ParameterBlock&
ParameterBlock::operator=(const ParameterBlock& other)
{
  if (this != &other)
  {
    type_ = other.type_;
    name_ = other.name_;
    if (other.value_ptr_)
      value_ptr_ = std::make_unique<Varying>(*other.value_ptr_);
    parameters_ = other.parameters_;
    error_origin_scope_ = other.error_origin_scope_;
  }

  return *this;
}

ParameterBlock::ParameterBlock(ParameterBlock&& other) noexcept
{
  std::swap(type_, other.type_);
  std::swap(name_, other.name_);
  std::swap(value_ptr_, other.value_ptr_);
  std::swap(parameters_, other.parameters_);
  std::swap(error_origin_scope_, other.error_origin_scope_);
}

ParameterBlock&
ParameterBlock::operator=(ParameterBlock&& other) noexcept
{
  if (this != &other)
  {
    std::swap(type_, other.type_);
    std::swap(name_, other.name_);
    std::swap(value_ptr_, other.value_ptr_);
    std::swap(parameters_, other.parameters_);
    std::swap(error_origin_scope_, other.error_origin_scope_);
  }

  return *this;
}

// Accessors
ParameterBlockType
ParameterBlock::GetType() const
{
  return type_;
}

bool
ParameterBlock::IsScalar() const
{
  return (type_ >= ParameterBlockType::BOOLEAN and type_ <= ParameterBlockType::INTEGER);
}
std::string
ParameterBlock::GetTypeName() const
{
  return ParameterBlockTypeName(type_);
}
std::string
ParameterBlock::GetName() const
{
  return name_;
}

const Varying&
ParameterBlock::GetValue() const
{
  switch (this->GetType())
  {
    case ParameterBlockType::BOOLEAN:
    case ParameterBlockType::FLOAT:
    case ParameterBlockType::STRING:
    case ParameterBlockType::INTEGER:
    case ParameterBlockType::USER_DATA:
    {
      if (value_ptr_ == nullptr)
        throw std::runtime_error(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) +
                                 ": Uninitialized Varying value for block " + this->GetName());
      return *value_ptr_;
    }
    default:
      throw std::logic_error(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) + ":\"" +
                             this->GetName() +
                             "\""
                             " Called for block of type " +
                             ParameterBlockTypeName(this->GetType()) + " which has no value.");
  }
}

size_t
ParameterBlock::GetNumParameters() const
{
  return parameters_.size();
}

const std::vector<ParameterBlock>&
ParameterBlock::GetParameters() const
{
  return parameters_;
}

bool
ParameterBlock::HasValue() const
{
  return value_ptr_ != nullptr;
}

// Mutators

void
ParameterBlock::ChangeToArray()
{
  const std::string fname = __PRETTY_FUNCTION__;
  if (parameters_.empty())
  {
    type_ = ParameterBlockType::ARRAY;
    return;
  }

  const auto& first_param = parameters_.front();
  for (const auto& param : parameters_)
    if (param.GetType() != first_param.GetType())
      throw std::logic_error(error_origin_scope_ + fname +
                             ": Cannot change ParameterBlock to "
                             "array. It has existing parameters and they are not of the same"
                             "type.");

  type_ = ParameterBlockType::ARRAY;
}

// NOLINTBEGIN(misc-no-recursion)
void
ParameterBlock::SetErrorOriginScope(const std::string& scope)
{
  error_origin_scope_ = scope;
  for (auto& param : parameters_)
    param.SetErrorOriginScope(scope);
}
// NOLINTEND(misc-no-recursion)

void
ParameterBlock::RequireBlockTypeIs(ParameterBlockType type) const
{
  if (GetType() != type)
    throw std::logic_error(error_origin_scope_ + ":" + GetName() + " Is required to be of type " +
                           ParameterBlockTypeName(type) + " but is " +
                           ParameterBlockTypeName(GetType()));
}

void
ParameterBlock::RequireParameter(const std::string& param_name) const
{
  if (not Has(param_name))
    throw std::logic_error(error_origin_scope_ + ":" + GetName() +
                           " Is required to have parameter " + param_name);
}

void
ParameterBlock::AddParameter(ParameterBlock block)
{
  for (const auto& param : parameters_)
    if (param.GetName() == block.GetName())
      throw std::invalid_argument(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) +
                                  ": Attempting to add duplicate parameter " + param.GetName() +
                                  " to "
                                  "block " +
                                  this->GetName());
  parameters_.push_back(std::move(block));

  SortParameters();
}

void
ParameterBlock::SortParameters()
{
  struct AlphabeticFunctor
  {
    bool operator()(const ParameterBlock& paramA, const ParameterBlock& paramB)
    {
      return paramA.GetName() < paramB.GetName();
    }
  };

  struct AlphabeticNumericFunctor
  {
    bool operator()(const ParameterBlock& paramA, const ParameterBlock& paramB)
    {
      return std::stoi(paramA.GetName()) < std::stoi(paramB.GetName());
    }
  };

  // The different functor here is necessary when the parameters are guaranteed
  // to have integer names that were converted to strings. It never showed up
  // because we were not testing with enough values. Essentially "11" < "2" in
  // the realm of strings but not in integer world.
  if (this->GetType() != ParameterBlockType::ARRAY)
    std::sort(parameters_.begin(), parameters_.end(), AlphabeticFunctor());
  else
    std::sort(parameters_.begin(), parameters_.end(), AlphabeticNumericFunctor());
}

bool
ParameterBlock::Has(const std::string& param_name) const
{
  return std::any_of(parameters_.begin(),
                     parameters_.end(),
                     [&param_name](const ParameterBlock& param)
                     { return param.name_ == param_name; });
}

ParameterBlock&
ParameterBlock::GetParam(const std::string& param_name)
{
  for (auto& param : parameters_)
    if (param.GetName() == param_name)
      return param;

  throw std::out_of_range(error_origin_scope_ + ":" + std::string(__PRETTY_FUNCTION__) +
                          ": Parameter \"" + param_name + "\" not present in block");
}

ParameterBlock&
ParameterBlock::GetParam(size_t index)
{
  try
  {
    return parameters_.at(index);
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) +
                            ": Parameter with index " + std::to_string(index) +
                            " not present in block");
  }
}

const ParameterBlock&
ParameterBlock::GetParam(const std::string& param_name) const
{
  for (const auto& param : parameters_)
    if (param.GetName() == param_name)
      return param;

  throw std::out_of_range(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) +
                          ": Parameter \"" + param_name + "\" not present in block");
}

const ParameterBlock&
ParameterBlock::GetParam(size_t index) const
{
  try
  {
    return parameters_.at(index);
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) +
                            ": Parameter with index " + std::to_string(index) +
                            " not present in block");
  }
}

//  NOLINTBEGIN(misc-no-recursion)
void
ParameterBlock::RecursiveDumpToString(std::string& outstr, const std::string& offset) const
{
  outstr += offset + this->GetName() + " = \n";
  outstr += offset + "{\n";

  if (HasValue())
    outstr += value_ptr_->PrintStr();

  for (const auto& param : parameters_)
  {

    switch (param.GetType())
    {
      case ParameterBlockType::BOOLEAN:
      {
        outstr += offset + "  " + param.GetName() + " = ";
        const bool value = param.GetValue().GetBoolValue();
        outstr += std::string(value ? "true" : "false") + ",\n";
        break;
      }
      case ParameterBlockType::FLOAT:
      {
        outstr += offset + "  " + param.GetName() + " = ";
        const double value = param.GetValue().GetFloatValue();
        outstr += std::to_string(value) + ",\n";
        break;
      }
      case ParameterBlockType::STRING:
      {
        outstr += offset + "  " + param.GetName() + " = ";
        const auto& value = param.GetValue().GetStringValue();
        outstr += "\"" + value + "\",\n";
        break;
      }
      case ParameterBlockType::INTEGER:
      {
        outstr += offset + "  " + param.GetName() + " = ";
        const int64_t value = param.GetValue().GetIntegerValue();
        outstr += std::to_string(value) + ",\n";
        break;
      }
      case ParameterBlockType::ARRAY:
      case ParameterBlockType::BLOCK:
      {
        param.RecursiveDumpToString(outstr, offset + "  ");
        break;
      }
      default:
        break;
    }
  } // for parameter

  outstr += offset + "}\n";
}
// NOLINTEND(misc-no-recursion)

//  NOLINTBEGIN(misc-no-recursion)
void
ParameterBlock::RecursiveDumpToJSON(std::string& outstr) const
{
  if (HasValue())
  {
    outstr += value_ptr_->PrintStr(false);
    return;
  }

  outstr += (this->GetType() == ParameterBlockType::ARRAY ? "[" : "{");
  for (const auto& param : parameters_)
  {

    switch (param.GetType())
    {
      case ParameterBlockType::BOOLEAN:
      {
        outstr += "\"" + param.GetName() + "\" = ";
        const bool value = param.GetValue().GetBoolValue();
        outstr += std::string(value ? "true" : "false") + ",\n";
        break;
      }
      case ParameterBlockType::FLOAT:
      {
        outstr += "\"" + param.GetName() + "\" = ";
        const double value = param.GetValue().GetFloatValue();
        outstr += std::to_string(value) + ",\n";
        break;
      }
      case ParameterBlockType::STRING:
      {
        outstr += "\"" + param.GetName() + "\" = ";
        const auto& value = param.GetValue().GetStringValue();
        outstr += "\"" + value + "\",\n";
        break;
      }
      case ParameterBlockType::INTEGER:
      {
        outstr += "\"" + param.GetName() + "\" = ";
        const int64_t value = param.GetValue().GetIntegerValue();
        outstr += std::to_string(value) + ",\n";
        break;
      }
      case ParameterBlockType::ARRAY:
      case ParameterBlockType::BLOCK:
      {
        param.RecursiveDumpToJSON(outstr);
        break;
      }
      default:
        break;
    }
  } // for parameter
  outstr += (this->GetType() == ParameterBlockType::ARRAY ? "]" : "}");
}
// NOLINTEND(misc-no-recursion)

} // namespace opensn
