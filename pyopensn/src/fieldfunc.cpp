// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "pyapi.hpp"

#include <memory>
#include <string>
#include <vector>

#include "framework/field_functions/field_function.h"
#include "framework/field_functions/field_function_grid_based.h"

namespace opensn {

// Wrap field functions
static void wrap_field_function(py::module & ffunc) {
    // field function
    auto field_function = py::class_<FieldFunction, std::shared_ptr<FieldFunction>>(
        ffunc,
        "FieldFunction",
        R"(
        Field function.

        Wrapper of :cpp:class:`opensn::FieldFunction`.
        )"
    );
    // field function grid based
    auto field_function_grid_based = py::class_<FieldFunctionGridBased, std::shared_ptr<FieldFunctionGridBased>, FieldFunction>(
        ffunc,
        "FieldFunctionGridBased",
        R"(
        Field function grid based.

        Wrapper of :cpp:class:`opensn::FieldFunctionGridBased`.
        )"
    );
    field_function_grid_based.def_static(
        "ExportMultipleToVTK",
        [](py::list ff_list, const std::string & base_name) {
            std::vector<std::shared_ptr<const FieldFunctionGridBased>> cpp_ff_list;
            cpp_ff_list.reserve(ff_list.size());
            for (py::handle item : ff_list) {
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
}

// Wrap the field function components of OpenSn
void py_ffunc(py::module & pyopensn) {
    py::module ffunc = pyopensn.def_submodule("fieldfunc", "Field function module.");
    wrap_field_function(ffunc);
}

}  // namespace opensn
