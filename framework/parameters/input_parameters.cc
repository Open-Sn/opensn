// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/parameters/input_parameters.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <sstream>
#include <algorithm>
#include <utility>

namespace opensn
{

namespace
{

std::string
InputErrorStr(const std::string& error_scope,
              const std::string& object_type,
              const std::string& err_msg)
{
  return (error_scope.empty() ? "" : error_scope + "\n") + "Input error: " + object_type + "\n" +
         err_msg;
}

std::string
ParamNotPresentErrorStr(const std::string& function_name, const std::string& param_name)
{
  return std::string(__PRETTY_FUNCTION__) + ": Parameter \"" + param_name +
         "\" not present in list of parameters.";
}

} // namespace

const std::vector<std::string> InputParameters::system_ignored_param_names_ = {"obj_type"};

InputParameters&
InputParameters::operator+=(InputParameters other)
{
  for (const auto& param : other)
    AddParameter(param);

  // Copy maps
  {
    auto& other_map = other.parameter_class_tags_;
    auto& this_map = parameter_class_tags_;
    for (const auto& [param_name, tag] : other_map)
    {
      OpenSnLogicalErrorIf(this_map.count(param_name) != 0, "Duplicate tags detected.");
      this_map[param_name] = tag;
    }
  }
  {
    auto& other_map = other.parameter_doc_string_;
    auto& this_map = parameter_doc_string_;
    for (const auto& [param_name, tag] : other_map)
    {
      OpenSnLogicalErrorIf(this_map.count(param_name) != 0, "Duplicate tags detected.");
      this_map[param_name] = tag;
    }
  }
  {
    auto& other_map = other.deprecation_warning_tags_;
    auto& this_map = deprecation_warning_tags_;
    for (const auto& [param_name, tag] : other_map)
    {
      OpenSnLogicalErrorIf(this_map.count(param_name) != 0, "Duplicate tags detected.");
      this_map[param_name] = tag;
    }
  }
  {
    auto& other_map = other.deprecation_error_tags_;
    auto& this_map = deprecation_error_tags_;
    for (const auto& [param_name, tag] : other_map)
    {
      OpenSnLogicalErrorIf(this_map.count(param_name) != 0, "Duplicate tags detected.");
      this_map[param_name] = tag;
    }
  }
  {
    auto& other_map = other.renamed_error_tags_;
    auto& this_map = renamed_error_tags_;
    for (const auto& [param_name, tag] : other_map)
    {
      OpenSnLogicalErrorIf(this_map.count(param_name) != 0, "Duplicate tags detected.");
      this_map[param_name] = tag;
    }
  }
  {
    auto& other_map = other.type_mismatch_allowed_tags_;
    auto& this_map = type_mismatch_allowed_tags_;
    for (const auto& [param_name, tag] : other_map)
    {
      OpenSnLogicalErrorIf(this_map.count(param_name) != 0, "Duplicate tags detected.");
      this_map[param_name] = tag;
    }
  }
  {
    auto& other_map = other.parameter_link_;
    auto& this_map = parameter_link_;
    for (const auto& [param_name, tag] : other_map)
    {
      OpenSnLogicalErrorIf(this_map.count(param_name) != 0, "Duplicate tags detected.");
      this_map[param_name] = tag;
    }
  }
  {
    auto& other_map = other.constraint_tags_;
    auto& this_map = constraint_tags_;
    for (auto& [param_name, tag] : other_map)
    {
      OpenSnLogicalErrorIf(this_map.count(param_name) != 0, "Duplicate tags detected.");
      this_map[param_name] = std::move(tag);
    }
  }

  return *this;
}

void
InputParameters::SetObjectType(const std::string& obj_type)
{
  class_name_ = obj_type;
}

std::string
InputParameters::GetObjectType() const
{
  return class_name_;
}

void
InputParameters::LinkParameterToBlock(const std::string& param_name, const std::string& block_name)
{
  OpenSnInvalidArgumentIf(not this->Has(param_name),
                          "Parameter \"" + param_name + "\" not present in block");
  parameter_link_[param_name] = block_name;
}

std::string
InputParameters::GetParameterDocumentationLink(const std::string& param_name) const
{
  if (parameter_link_.count(param_name) > 0)
    return parameter_link_.at(param_name);
  return {};
}

std::string
InputParameters::GetParameterDocString(const std::string& param_name)
{
  OpenSnInvalidArgumentIf(parameter_doc_string_.count(param_name) == 0,
                          "Invalid parameter \"" + param_name + "\".");
  return parameter_doc_string_.at(param_name);
}

bool
InputParameters::IsParameterIgnored(const std::string& param_name)
{
  bool ignored = false;

  {
    auto& list = system_ignored_param_names_;
    if (std::find(list.begin(), list.end(), param_name) != list.end())
      ignored = true;
  }

  return ignored;
}

void
InputParameters::AddOptionalParameterBlock(const std::string& name,
                                           const ParameterBlock& block,
                                           const std::string& doc_string)
{
  auto new_block = block;
  new_block.SetBlockName(name);
  AddParameter(new_block);
  parameter_class_tags_[name] = InputParameterTag::OPTIONAL;
  parameter_doc_string_[name] = doc_string;
}

void
InputParameters::AddOptionalParameterArray(const std::string& name,
                                           const std::vector<ParameterBlock>& array,
                                           const std::string& doc_string)
{
  ParameterBlock new_block(name);
  new_block.ChangeToArray();
  for (auto& block : array)
    new_block.AddParameter(block);

  AddParameter(new_block);
  parameter_class_tags_[name] = InputParameterTag::OPTIONAL;
  parameter_doc_string_[name] = doc_string;
}

void
InputParameters::AddRequiredParameterBlock(const std::string& name, const std::string& doc_string)
{
  ParameterBlock new_block(name);
  AddParameter(new_block);
  parameter_class_tags_[name] = InputParameterTag::REQUIRED;
  parameter_doc_string_[name] = doc_string;
}

void
InputParameters::AddRequiredParameterArray(const std::string& name, const std::string& doc_string)
{
  ParameterBlock new_block(name);
  new_block.ChangeToArray();
  AddParameter(new_block);
  parameter_class_tags_[name] = InputParameterTag::REQUIRED;
  parameter_doc_string_[name] = doc_string;
}

void
InputParameters::AssignParameters(const ParameterBlock& params)
{
  param_block_at_assignment_ = params;
  std::stringstream err_stream;

  if (log.GetVerbosity() >= 2)
    log.Log0Verbose2() << "Number of parameters " << params.GetNumParameters();

  // Check required parameters
  // Loops over all input-parameters that have been
  // classed as being required. Input-parameters that
  // have any form of deprecation is ignored.
  for (const auto& [param_index, tag] : parameter_class_tags_)
  {
    if (tag != InputParameterTag::REQUIRED)
      continue;

    const auto& req_param = GetParam(param_index);
    const auto& req_param_name = req_param.GetName();

    if (deprecation_warning_tags_.count(req_param_name) > 0 or
        deprecation_error_tags_.count(req_param_name) > 0 or
        renamed_error_tags_.count(req_param_name) > 0)
      continue;

    if (not params.Has(req_param_name))
      err_stream << "Required param \"" << req_param_name
                 << "\" not supplied.\ndoc-string: " << GetParameterDocString(req_param_name)
                 << "\nEnsure the parameter given is supplied or not nil";
  }

  if (not err_stream.str().empty())
    throw std::invalid_argument(
      InputErrorStr(GetErrorOriginScope(), GetObjectType(), err_stream.str()));

  // Check unused parameters
  // Loops over all candidate-parameters and
  // checks whether they have an assignable
  // input-parameter or if they have been renamed.
  {
    for (const auto& param : params.GetParameters())
    {
      const auto& param_name = param.GetName();
      if (IsParameterIgnored(param_name))
        continue;
      if (not this->Has(param_name))
        err_stream << "Invalid param \"" << param_name << "\" supplied.\n";
      else if (renamed_error_tags_.count(param_name) > 0)
      {
        err_stream << "Invalid param \"" << param_name << "\" supplied. ";
        err_stream << "The parameter has been renamed. ";
        err_stream << renamed_error_tags_.at(param_name);
      }
    }

    if (not err_stream.str().empty())
      throw std::invalid_argument(
        InputErrorStr(GetErrorOriginScope(), GetObjectType(), err_stream.str()));
  }

  // Check deprecation warnings
  // Loops over all candidate-parameters and
  // checks whether they have deprecation warnings.
  {
    const auto& dep_warns = deprecation_warning_tags_;
    for (const auto& param : params.GetParameters())
    {
      const auto& param_name = param.GetName();

      if (IsParameterIgnored(param_name))
        continue;

      if (this->Has(param_name) and (dep_warns.count(param_name) > 0))
        log.Log0Warning() << "Parameter \"" << param_name << "\" has been deprecated "
                          << "and will be removed soon.\n"
                          << dep_warns.at(param_name);
    }
  }

  // Check deprecation errors
  // Loops over all candidate-parameters and
  // checks whether they have deprecation errors.
  {
    const auto& dep_errs = deprecation_error_tags_;
    for (const auto& param : params.GetParameters())
    {
      const auto& param_name = param.GetName();

      if (IsParameterIgnored(param_name))
        continue;

      if (this->Has(param_name) and (dep_errs.count(param_name) > 0))
      {
        log.Log0Error() << "Parameter \"" << param_name << "\" has been deprecated.\n"
                        << dep_errs.at(param_name);
        Exit(EXIT_FAILURE);
      }
    }
  }

  // Now attempt to assign values
  for (auto& param : params.GetParameters())
  {
    const auto& param_name = param.GetName();

    if (IsParameterIgnored(param_name))
      continue;

    auto& input_param = GetParam(param_name);

    // Check types match
    if (param.GetType() != input_param.GetType())
    {
      if (type_mismatch_allowed_tags_.count(param_name) == 0)
      {
        err_stream << "Invalid parameter type \"" << ParameterBlockTypeName(param.GetType())
                   << "\" for parameter \"" << param_name << "\". Expecting type \""
                   << ParameterBlockTypeName(input_param.GetType()) << "\".\n"
                   << "doc-string: " << GetParameterDocString(param_name);
        continue;
      } // if not mismatch allowed
    }   // if type mismatch

    // Check constraint
    if (constraint_tags_.count(input_param.GetName()) != 0)
    {
      const auto& constraint = constraint_tags_.at(input_param.GetName());
      if (not constraint->IsAllowable(param.GetValue()))
      {
        err_stream << constraint->OutOfRangeString(input_param.GetName(), param.GetValue());
        err_stream << "\n";
        continue;
      }
    } // if constraint

    if (log.GetVerbosity() >= 2)
      log.Log0Verbose2() << "Setting parameter " << param_name;
    input_param = param;
    parameter_valid_[param_name] = true;
  } // for input params

  if (not err_stream.str().empty())
    throw std::invalid_argument(
      InputErrorStr(GetErrorOriginScope(), GetObjectType(), err_stream.str()));
}

void
InputParameters::MarkParamaterDeprecatedWarning(const std::string& param_name,
                                                const std::string& deprecation_message)
{
  if (Has(param_name))
    deprecation_warning_tags_[param_name] = deprecation_message;
  else
    throw std::logic_error(ParamNotPresentErrorStr(__PRETTY_FUNCTION__, param_name));
}

void
InputParameters::MarkParamaterDeprecatedError(const std::string& param_name,
                                              const std::string& deprecation_message)
{
  if (Has(param_name))
    deprecation_error_tags_[param_name] = deprecation_message;
  else
    throw std::logic_error(ParamNotPresentErrorStr(__PRETTY_FUNCTION__, param_name));
}

void
InputParameters::MarkParamaterRenamed(const std::string& param_name,
                                      const std::string& renaming_description)
{
  if (Has(param_name))
    renamed_error_tags_[param_name] = renaming_description;
  else
    throw std::logic_error(ParamNotPresentErrorStr(__PRETTY_FUNCTION__, param_name));
}

void
InputParameters::ConstrainParameterRange(const std::string& param_name,
                                         std::shared_ptr<AllowableRange> allowable_range)
{
  if (Has(param_name))
  {
    const auto& param_type = GetParam(param_name).GetType();
    OpenSnInvalidArgumentIf(
      param_type == ParameterBlockType::BLOCK or param_type == ParameterBlockType::ARRAY,
      std::string("Parameter \"") + param_name + "\" is of type " +
        ParameterBlockTypeName(param_type) + " to which constraints cannot be applied");
    constraint_tags_[param_name] = allowable_range;
  }
  else
    throw std::logic_error(ParamNotPresentErrorStr(__PRETTY_FUNCTION__, param_name));
}

void
InputParameters::SetParameterTypeMismatchAllowed(const std::string& param_name)
{
  OpenSnInvalidArgumentIf(not Has(param_name), "Parameter \"" + param_name + "\" not present.");
  type_mismatch_allowed_tags_[param_name] = true;
}

void
InputParameters::DumpParameters() const
{
  log.Log() << "CLASS_NAME " << class_name_;

  log.Log() << "DESCRIPTION_BEGIN";
  std::cout << GetGeneralDescription() << "\n";
  log.Log() << "DESCRIPTION_END\n";

  log.Log() << "DOC_GROUP " << doc_group_;
  const std::string sp2 = "  ";
  const std::string sp4 = "    ";
  const auto params = GetParameters();
  for (const auto& param : params)
  {
    const auto& param_name = param.GetName();
    log.Log() << sp2 << "PARAM_BEGIN " << param_name;

    const auto type = param.GetType();

    log.Log() << sp4 << "TYPE " << ParameterBlockTypeName(type);

    if (parameter_class_tags_.at(param_name) == InputParameterTag::OPTIONAL)
    {
      log.Log() << sp4 << "TAG OPTIONAL";
      if (type != ParameterBlockType::BLOCK and type != ParameterBlockType::ARRAY)
        log.Log() << sp4 << "DEFAULT_VALUE " << param.GetValue().PrintStr();
      else if (type == ParameterBlockType::ARRAY)
      {
        std::stringstream outstr;
        outstr << sp4 << "DEFAULT_VALUE ";
        for (size_t k = 0; k < param.GetNumParameters(); ++k)
        {
          const auto& sub_param = param.GetParam(k);
          outstr << sub_param.GetValue().PrintStr() << ", ";
        }
        log.Log() << outstr.str();
      }
    }
    else
      log.Log() << sp4 << "TAG REQUIRED";

    if (constraint_tags_.count(param_name) != 0)
      log.Log() << sp4 << "CONSTRAINTS " << constraint_tags_.at(param_name)->PrintRange();

    if (parameter_doc_string_.count(param_name) != 0)
    {
      log.Log() << sp4 << "DOC_STRING_BEGIN";
      std::cout << parameter_doc_string_.at(param_name) << "\n";
      log.Log() << sp4 << "DOC_STRING_END";
    }

    const auto& linkage = GetParameterDocumentationLink(param_name);
    if (not linkage.empty())
    {
      log.Log() << sp4 << "LINKS " << linkage;
    }

    log.Log() << sp2 << "PARAM_END";
  }
}

bool
InputParameters::IsParameterValid(const std::string& param_name) const
{
  auto it = parameter_valid_.find(param_name);
  if (it != parameter_valid_.end())
    return it->second;
  else
    return false;
}

} // namespace opensn
