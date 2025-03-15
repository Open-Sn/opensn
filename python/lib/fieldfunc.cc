// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "python/lib/functor.h" // temporary, see the included header for more details!
#include "framework/field_functions/field_function.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/field_functions/interpolation/ffinter_point.h"
#include "framework/field_functions/interpolation/ffinter_line.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace opensn
{

// dictionary to field function interpolation operation type
static std::map<std::string, FieldFunctionInterpolationOperation> ff_op_type_map{
  {"sum", FieldFunctionInterpolationOperation::OP_SUM},
  {"avg", FieldFunctionInterpolationOperation::OP_AVG},
  {"max", FieldFunctionInterpolationOperation::OP_MAX},
  {"sum_func", FieldFunctionInterpolationOperation::OP_SUM_FUNC},
  {"avg_func", FieldFunctionInterpolationOperation::OP_AVG_FUNC},
  {"max_func", FieldFunctionInterpolationOperation::OP_MAX_FUNC}};

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
    "ExportMultipleToVTK",
    [](py::list& ff_list, const std::string& base_name)
    {
      std::vector<std::shared_ptr<const FieldFunctionGridBased>> cpp_ff_list;
      cpp_ff_list.reserve(ff_list.size());
      for (py::handle item : ff_list)
      {
        cpp_ff_list.push_back(item.cast<std::shared_ptr<const FieldFunctionGridBased>>());
      }
      FieldFunctionGridBased::ExportMultipleToVTK(base_name, cpp_ff_list);
    },
    R"(
      Export a list of "field function grid based" to VTK format.

      Parameters
      ----------
      ff_list: List[opensn.FieldFunctionGridBased]
          List of "field function grid based" to export.
      base_name: str
          Base name.
    )",
    py::arg("ff_list"), py::arg("base_name")
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

    Wrapper of :cpp:class:`opensn::FieldFunctionInterpolation`.
    )"
  );
  field_func_interp.def(
    "AddFieldFunction",
    &FieldFunctionInterpolation::AddFieldFunction,
    R"(
    Add a field function to the list.
    )",
    py::arg("ff")
  );
  field_func_interp.def(
    "Initialize",
    &FieldFunctionInterpolation::Initialize,
    R"(
    ???
    )"
  );
  field_func_interp.def(
    "Execute",
    &FieldFunctionInterpolation::Execute,
    R"(
    ???
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
  field_func_interp.def_static(
    "GetFieldFunctionByName",
    [](const std::string& ff_name)
    {
      // get list of suitable field functions
      py::list matched_ff;
      for (std::shared_ptr<FieldFunction>& ff_ptr : field_function_stack)
      {
        if (ff_ptr->GetName() == ff_name)
        {
          matched_ff.append(ff_ptr);
        }
      }
      return matched_ff;
    },
    R"(
    Get the list of field functions matching a given name.

    This function returns a list of field functions whose names match the given argument. The list
    may be empty or contain multiple elements.

    Parameters
    ----------
    ff_name: str
        Field function name
    )",
    py::arg("ff_name")
  );

  // field function interpolation point
  auto field_func_interp_point = py::class_<FieldFunctionInterpolationPoint,
                                            std::shared_ptr<FieldFunctionInterpolationPoint>,
                                            FieldFunctionInterpolation>(
    ffunc,
    "FieldFunctionInterpolationPoint",
    R"(
    Line based interpolation function.
    ??? (same docuumentation as FieldFunctionInterpolationLine)

    Wrapper of :cpp:class:`opensn::FieldFunctionInterpolationPoint`.
    )"
  );
  field_func_interp_point.def(
    py::init(
      [](void)
      {
        return FieldFunctionInterpolationPoint::Create();
      }
    ),
    "Default constructor."
  );
  field_func_interp_point.def(
    "GetPointValue",
    &FieldFunctionInterpolationPoint::GetPointValue,
    R"(
    ???
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

    Wrapper of :cpp:class:`opensn::FieldFunctionInterpolationLine`.
    )"
  );
  field_func_interp_line.def(
    py::init(
      [](void)
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
    point: List[float]
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
    point: List[float]
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
    A line based interpolation function.

    Wrapper of :cpp:class:`opensn::FieldFunctionInterpolationVolume`.
    )"
  );
  field_func_interp_volume.def(
    py::init(
      [](void)
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
      self.SetOperationType(ff_op_type_map.at(op_type));
    },
    R"(
    Set operation type.

    Parameters
    ----------
    op_type: {'sum', 'avg', 'max', 'sum_func', 'avg_func', 'max_func'}
        Operation type.
    )",
    py::arg("op_type")
  );
  field_func_interp_volume.def(
    "SetOperationFunction",
    [](FieldFunctionInterpolationVolume& self, std::shared_ptr<PySMFunction> function)
    {
      self.SetOperationFunction(function);
    },
    R"(
    ???

    Parameters
    ----------
    function: pyopensn.math.ScalarMaterialFunction
        ???
    )",
    py::arg("function")
  );
  field_func_interp_volume.def(
    "GetValue",
    &FieldFunctionInterpolationVolume::GetValue,
    R"(
    ???
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
