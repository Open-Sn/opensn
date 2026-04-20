// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/field_functions/field_function.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/field_functions/interpolation/ffinter_point.h"
#include "framework/field_functions/interpolation/ffinter_line.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include <pybind11/functional.h>
#include <algorithm>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace opensn
{

namespace
{

// dictionary to field function interpolation operation type
std::map<std::string, FieldFunctionInterpolationOperation> ff_op_type_map{
  {"sum", FieldFunctionInterpolationOperation::OP_SUM},
  {"avg", FieldFunctionInterpolationOperation::OP_AVG},
  {"max", FieldFunctionInterpolationOperation::OP_MAX},
  {"min", FieldFunctionInterpolationOperation::OP_MIN},
  {"sum_func", FieldFunctionInterpolationOperation::OP_SUM_FUNC},
  {"avg_func", FieldFunctionInterpolationOperation::OP_AVG_FUNC},
  {"max_func", FieldFunctionInterpolationOperation::OP_MAX_FUNC},
  {"min_func", FieldFunctionInterpolationOperation::OP_MIN_FUNC}};

std::string
OpToString()
{
  std::vector<std::string> keys;
  keys.reserve(ff_op_type_map.size());
  for (const auto& [key, _] : ff_op_type_map)
    keys.push_back(key);
  std::sort(keys.begin(), keys.end());
  std::ostringstream oss;
  for (size_t i = 0; i < keys.size(); ++i)
  {
    if (i > 0)
      oss << ", ";
    oss << keys[i];
  }
  return oss.str();
}

} // namespace

// Wrap field functions
void
WrapFieldFunction(py::module& ffunc)
{
  // clang-format off
  // field function
  auto field_func = py::class_<FieldFunction, std::shared_ptr<FieldFunction>>(
    ffunc,
    "FieldFunction",
    R"(
    Field function.

    Wrapper of :cpp:class:`opensn::FieldFunction`.
    )"
  );
  // clang-format on
}

// Wrap field functions grid based
void
WrapFieldFunctionGridBased(py::module& ffunc)
{
  // clang-format off
  // field function grid based
  auto field_func_grid_based = py::class_<FieldFunctionGridBased,
                                          std::shared_ptr<FieldFunctionGridBased>,
                                          FieldFunction>(
    ffunc,
    "FieldFunctionGridBased",
    R"(
    Field function grid based.

    Wrapper of :cpp:class:`opensn::FieldFunctionGridBased`.
    )"
  );
  field_func_grid_based.def_static(
    "ExportMultipleToPVTU",
    [](py::list& ff_list, const std::string& base_name)
    {
      std::vector<std::shared_ptr<const FieldFunctionGridBased>> cpp_ff_list;
      cpp_ff_list.reserve(ff_list.size());
      for (py::handle item : ff_list)
      {
        cpp_ff_list.push_back(item.cast<std::shared_ptr<const FieldFunctionGridBased>>());
      }
      FieldFunctionGridBased::ExportMultipleToPVTU(base_name, cpp_ff_list);
    },
    R"(
      Export a list of "field function grid based" to parallel VTU format.

      Parameters
      ----------
      ff_list: List[opensn.FieldFunctionGridBased]
          List of "field function grid based" to export.
      base_name: str
          Base name.
    )",
    py::arg("ff_list"), py::arg("base_name")
  );
  field_func_grid_based.def(
    "CanUpdate",
    &FieldFunctionGridBased::CanUpdate,
    R"(
    Return whether this field function can currently refresh itself from its owning problem.

    This returns ``False`` if the field function has no update callback or if the owning problem
    has already been destroyed.
    )"
  );
  field_func_grid_based.def(
    "Update",
    &FieldFunctionGridBased::Update,
    R"(
    Refresh this field function from its owning problem.

    Raises an error if the field function is not refreshable or if its owning problem has already
    been destroyed.
    )"
  );
  // clang-format on
}

// Wrap field function interpolation
void
WrapFieldFunctionInterpolation(py::module& ffunc)
{
  // clang-format off
  // field function interpolation
  auto field_func_interp = py::class_<FieldFunctionInterpolation,
                                      std::shared_ptr<FieldFunctionInterpolation>>(
    ffunc,
    "FieldFunctionInterpolation",
    R"(
    Base class for field-function interpolation objects.

    Interpolators are configured by assigning a field function and any
    geometry- or operation-specific settings, then calling ``Execute()``.
    ``Execute()`` rebuilds any required internal state from the current
    configuration each time it is called.

    Wrapper of :cpp:class:`opensn::FieldFunctionInterpolation`.
    )"
  );
  field_func_interp.def(
    "AddFieldFunction",
    &FieldFunctionInterpolation::AddFieldFunction,
    R"(
    Add a field function to this interpolator.

    Current point/line/volume interpolators support exactly one field function.
    This method only succeeds when no field function is currently assigned.
    It raises an error if one is already present. Use ``SetFieldFunction`` to
    replace the current one or ``ClearFieldFunctions`` to remove it first.

    Changing the assigned field function affects subsequent ``Execute()``
    calls; no separate initialization step is required.
    )",
    py::arg("ff")
  );
  field_func_interp.def(
    "SetFieldFunction",
    &FieldFunctionInterpolation::SetFieldFunction,
    R"(
    Replace the current field function with ``ff``.

    Unlike ``AddFieldFunction``, this method does not require the interpolator
    to be empty first.

    The next call to ``Execute()`` uses this field function.
    )",
    py::arg("ff")
  );
  field_func_interp.def(
    "ClearFieldFunctions",
    &FieldFunctionInterpolation::ClearFieldFunctions,
    R"(
    Remove all field functions from this interpolator.

    ``Execute()`` will raise an error until a new field function is assigned.
    )"
  );
  field_func_interp.def(
    "Execute",
    &FieldFunctionInterpolation::Execute,
    R"(
    Execute the field function interpolator using the current configuration.

    This method rebuilds any required internal interpolation state each time it
    is called. No separate ``Initialize()`` step is required.
    )"
  );
  field_func_interp.def(
    "ExportToCSV",
    &FieldFunctionInterpolation::ExportToCSV,
    R"(
    Export field function interpolation to CSV files.

    Parameters
    ----------
    base_name: str
        Base name of the exported CSVs.
    )",
    py::arg("base_name")
  );
  // field function interpolation point
  auto field_func_interp_point = py::class_<FieldFunctionInterpolationPoint,
                                            std::shared_ptr<FieldFunctionInterpolationPoint>,
                                            FieldFunctionInterpolation>(
    ffunc,
    "FieldFunctionInterpolationPoint",
    R"(
    Interpolate the field function at a point.

    Configure the point with ``SetPointOfInterest(...)``, assign a field
    function, then call ``Execute()`` before reading the result with
    ``GetPointValue()``.

    Wrapper of :cpp:class:`opensn::FieldFunctionInterpolationPoint`.
    )"
  );
  field_func_interp_point.def(
    py::init(
      []()
      {
        return FieldFunctionInterpolationPoint::Create();
      }
    ),
    "Default constructor."
  );
  field_func_interp_point.def(
    "SetPointOfInterest",
    &FieldFunctionInterpolationPoint::SetPointOfInterest,
    R"(
    Set the point at which the field function will be evaluated.

    Parameters
    ----------
    point: pyopensn.math.Vector3
        Coordinates of the point of interest.
    )",
    py::arg("point")
  );
  field_func_interp_point.def(
    "GetPointValue",
    &FieldFunctionInterpolationPoint::GetPointValue,
    R"(
    Get the most recently computed point value.

    Call ``Execute()`` after changing the point or field function to refresh
    the stored result.
    )"
  );

  // field function interpolation line
  auto field_func_interp_line = py::class_<FieldFunctionInterpolationLine,
                                           std::shared_ptr<FieldFunctionInterpolationLine>,
                                           FieldFunctionInterpolation>(
    ffunc,
    "FieldFunctionInterpolationLine",
    R"(
    Line based interpolation function.

    Configure the line segment and number of points, assign a field function,
    then call ``Execute()``. ``Execute()`` rebuilds the line sampling data from
    the current configuration before evaluating the field.

    Wrapper of :cpp:class:`opensn::FieldFunctionInterpolationLine`.
    )"
  );
  field_func_interp_line.def(
    py::init(
      []()
      {
        return FieldFunctionInterpolationLine::Create();
      }
    ),
    "Default constructor."
  );
  field_func_interp_line.def(
    "SetInitialPoint",
    &FieldFunctionInterpolationLine::SetInitialPoint,
    R"(
    Set initial point.

    Parameters
    ----------
    point: pyopensn.math.Vector3
        Coordinates of the initial point.
    )",
    py::arg("point")
  );
  field_func_interp_line.def(
    "SetFinalPoint",
    &FieldFunctionInterpolationLine::SetFinalPoint,
    R"(
    Set final point.

    Parameters
    ----------
    point: pyopensn.math.Vector3
        Coordinates of the final point.
    )",
    py::arg("point")
  );
  field_func_interp_line.def(
    "SetNumberOfPoints",
    &FieldFunctionInterpolationLine::SetNumberOfPoints,
    R"(
    Set number of points.

    Parameters
    ----------
    number: int
        Number of points.
    )",
    py::arg("number")
  );

  // field function interpolation volume
  auto field_func_interp_volume = py::class_<FieldFunctionInterpolationVolume,
                                             std::shared_ptr<FieldFunctionInterpolationVolume>,
                                             FieldFunctionInterpolation>(
    ffunc,
    "FieldFunctionInterpolationVolume",
    R"(
    Volume based interpolation function.

    Configure the logical volume and operation, assign a field function, then
    call ``Execute()``. ``Execute()`` rebuilds the cell set from the current
    configuration before evaluating the field.

    Wrapper of :cpp:class:`opensn::FieldFunctionInterpolationVolume`.
    )"
  );
  field_func_interp_volume.def(
    py::init(
      []()
      {
        return FieldFunctionInterpolationVolume::Create();
      }
    ),
    "Default constructor."
  );
  field_func_interp_volume.def(
    "SetLogicalVolume",
    &FieldFunctionInterpolationVolume::SetLogicalVolume,
    "Set logical volume.",
    py::arg("lv")
  );
  field_func_interp_volume.def(
    "SetOperationType",
    [](FieldFunctionInterpolationVolume& self, const std::string& op_type)
    {
      const auto it = ff_op_type_map.find(op_type);
      if (it == ff_op_type_map.end())
        throw std::runtime_error("Unknown field-function operation '" + op_type +
                                 "'. Valid options are: " + OpToString());
      self.SetOperationType(it->second);
    },
    R"(
    Set operation type.

    Parameters
    ----------
    op_type: {'sum', 'avg', 'max', 'min', 'sum_func', 'avg_func', 'max_func', 'min_func'}
        Operation type.
    )",
    py::arg("op_type")
  );
  field_func_interp_volume.def(
    "SetOperationFunction",
    [](FieldFunctionInterpolationVolume& self, const ScalarMaterialFunction& function)
    {
      self.SetOperationFunction(function);
    },
    R"(
    Set the field function operation type to a custom scalar material function.

    Parameters
    ----------
    function: Callable[[float, int], float]
        A scalar material function that takes the field function value (float) and the
        block id (int) as parameters and returns a double.
    )",
    py::arg("function")
  );
  field_func_interp_volume.def(
    "GetValue",
    &FieldFunctionInterpolationVolume::GetValue,
    R"(
    Return the most recently computed interpolation value.

    Call ``Execute()`` after changing the logical volume, operation, operation
    function, or field function to refresh the stored result.
    )"
  );
  // clang-format on
}

// Wrap the field function components of OpenSn
void
py_ffunc(py::module& pyopensn)
{
  py::module ffunc = pyopensn.def_submodule("fieldfunc", "Field function module.");
  WrapFieldFunction(ffunc);
  WrapFieldFunctionGridBased(ffunc);
  WrapFieldFunctionInterpolation(ffunc);
}

} // namespace opensn
