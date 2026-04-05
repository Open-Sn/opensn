// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <cstdint>

namespace opensn
{
class FieldFunction;
class FieldFunctionGridBased;

enum class FieldFunctionInterpolationType : int
{
  SLICE = 1,
  LINE = 2,
  VOLUME = 3,
  POINT = 4
};

enum class FieldFunctionInterpolationOperation : int
{
  OP_SUM = 10,
  OP_AVG = 11,
  OP_MAX = 12,
  OP_MIN = 13,
  OP_SUM_FUNC = 14,
  OP_AVG_FUNC = 15,
  OP_MAX_FUNC = 16,
  OP_MIN_FUNC = 17,
};

enum class FieldFunctionInterpolationProperty : int
{
  PROBEPOINT = 0,
  SLICEPOINT = 1,
  SLICENORMAL = 2,
  SLICETANGENT = 3,
  SLICEBINORM = 4,
  OPERATION = 5,
  LOGICAL_VOLUME = 8,
  ADD_FIELD_FUNCTION = 9,
  SET_FIELD_FUNCTIONS = 10,
  FIRSTPOINT = 11,
  SECONDPOINT = 12,
  NUMBEROFPOINTS = 13,
  CUSTOM_ARRAY = 14,
};

/// Base class for field-function interpolation objects.
class FieldFunctionInterpolation
{
public:
  explicit FieldFunctionInterpolation(FieldFunctionInterpolationType type)
    : type_(type), ref_component_(0)
  {
  }

  virtual ~FieldFunctionInterpolation() = default;

  void ClearFieldFunctions() { field_function_.reset(); }

  void SetFieldFunction(std::shared_ptr<FieldFunction> ff)
  {
    ClearFieldFunctions();
    AddFieldFunction(std::move(ff));
  }

  void AddFieldFunction(std::shared_ptr<FieldFunction> ff)
  {
    auto gbff = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff);
    if (gbff)
    {
      if (field_function_)
        throw std::runtime_error(
          "FieldFunctionInterpolation currently supports exactly one field function. "
          "Use SetFieldFunction to replace it or ClearFieldFunctions before adding another.");
      field_function_ = std::move(gbff);
    }
    else
      throw std::runtime_error("Expected FieldFunctionGridBased field function");
  }

  FieldFunctionInterpolationType Type() const { return type_; }

  virtual void Execute() {};

  virtual void ExportToCSV(std::string base_name) const {};

  virtual void ExportToPython(std::string base_name) {};

protected:
  FieldFunctionInterpolationType type_;
  unsigned int ref_component_;
  std::shared_ptr<FieldFunctionGridBased> field_function_;
};

} // namespace opensn
