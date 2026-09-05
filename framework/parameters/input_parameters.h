// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/parameter_block.h"
#include "framework/data_types/allowable_range.h"
#include <map>
#include <memory>
#include <cstdint>

namespace opensn
{

enum class InputParameterTag
{
  NONE = 0,
  OPTIONAL = 1,
  REQUIRED = 2
};

/// Class for handling input parameters.
class InputParameters : public ParameterBlock
{
public:
  InputParameters() = default;
  InputParameters& operator+=(InputParameters other);

  /// Sets the object type string for more descriptive error messages.
  void SetObjectType(const std::string& obj_type);

  /// Returns the object type string.
  std::string GetObjectType() const;

  /// Returns the parameter's doc string.
  std::string GetParameterDocString(const std::string& param_name);

  template <typename T>
  void AddOptionalParameter(const std::string& name, T value, const std::string& doc_string)
  {
    AddParameter(name, value);
    parameter_class_tags_[name] = InputParameterTag::OPTIONAL;
    parameter_doc_string_[name] = doc_string;
  }

  /// Specialization for block type parameters.
  void AddOptionalParameterBlock(const std::string& name,
                                 const ParameterBlock& block,
                                 const std::string& doc_string);

  template <typename T>
  void AddOptionalParameterArray(const std::string& name,
                                 const std::vector<T>& array,
                                 const std::string& doc_string)
  {
    AddParameter(name, array);
    parameter_class_tags_[name] = InputParameterTag::OPTIONAL;
    parameter_doc_string_[name] = doc_string;
  }

  /// Specialization for block type parameters.
  void AddOptionalParameterArray(const std::string& name,
                                 const std::vector<ParameterBlock>& array,
                                 const std::string& doc_string);

  template <typename T>
  void AddRequiredParameter(const std::string& name, const std::string& doc_string)
  {
    AddParameter(name, Varying::DefaultValue<T>());
    parameter_class_tags_[name] = InputParameterTag::REQUIRED;
    parameter_doc_string_[name] = doc_string;
  }

  /// Specialization for block type parameters.
  void AddRequiredParameterBlock(const std::string& name, const std::string& doc_string);

  /// Specialization for array type parameters.
  void AddRequiredParameterArray(const std::string& name, const std::string& doc_string);

  template <typename T>
  void ChangeExistingParamToOptional(const std::string& name,
                                     T value,
                                     const std::string& doc_string = "")
  {
    auto& param = GetParam(name);
    param = ParameterBlock(name, value);
    parameter_class_tags_[name] = InputParameterTag::OPTIONAL;
    if (not doc_string.empty())
      parameter_doc_string_[name] = doc_string;
  }

  template <typename T>
  void ChangeExistingParamToRequired(const std::string& name, const std::string& doc_string = "")
  {
    auto& param = GetParam(name);
    param = ParameterBlock(name, Varying::DefaultValue<T>());
    parameter_class_tags_[name] = InputParameterTag::REQUIRED;
    if (not doc_string.empty())
      parameter_doc_string_[name] = doc_string;
  }

  /// Assigns parameters with thorough type checks, deprecation checks, unused parameter checks.
  void AssignParameters(const ParameterBlock& params);

  /**
   * Returns the raw parameter block used at assignment. This can be used to see if a user supplied
   * an optional parameter or not.
   */
  const ParameterBlock& GetParametersAtAssignment() const { return param_block_at_assignment_; }

  bool IsParameterValid(const std::string& param_name) const;

  /// Marks a parameters as deprecated but will only produce a warning.
  void MarkParameterDeprecatedWarning(const std::string& param_name,
                                      const std::string& deprecation_message = "");

  /// Marks a parameters as deprecated and will produce an error if the parameter is specified.
  void MarkParameterDeprecatedError(const std::string& param_name,
                                    const std::string& deprecation_message = "");

  /// Marks a parameters as renamed and will produce an error if the parameter is specified.
  void MarkParameterRenamed(const std::string& param_name, const std::string& renaming_description);

  /// Creates a range based constraint for a given parameter.
  void ConstrainParameterRange(const std::string& param_name,
                               std::shared_ptr<AllowableRange> allowable_range);

  /// Sets a tag for the given parameter that will allow its type to be mismatched upon assignment.
  void SetParameterTypeMismatchAllowed(const std::string& param_name);

private:
  /// Backing storage for the object type string set via SetObjectType().
  std::string class_name_;
  std::map<std::string, InputParameterTag> parameter_class_tags_;
  std::map<std::string, std::string> parameter_doc_string_;
  std::map<std::string, bool> parameter_valid_;
  std::map<std::string, std::string> deprecation_warning_tags_;
  std::map<std::string, std::string> deprecation_error_tags_;
  std::map<std::string, std::string> renamed_error_tags_;
  std::map<std::string, bool> type_mismatch_allowed_tags_;

  std::map<std::string, std::shared_ptr<AllowableRange>> constraint_tags_;

  ParameterBlock param_block_at_assignment_;

private:
  using ParameterBlock::AddParameter;
};

/**
 * Constructs a `T` from user-supplied parameters, using `T::GetInputParameters()` as the
 * parameter input. `obj_type` is used only to make error messages more descriptive.
 */
template <typename T>
std::shared_ptr<T>
CreateObject(const std::string& obj_type, const ParameterBlock& params)
{
  auto input_params = T::GetInputParameters();
  input_params.SetObjectType(obj_type);
  input_params.SetErrorOriginScope(obj_type);
  input_params.AssignParameters(params);
  return std::make_shared<T>(input_params);
}

} // namespace opensn
