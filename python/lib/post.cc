// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/cross_section_sensitivity_postprocessor.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/volume_postprocessor.h"
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace opensn
{

// Wrap post processors
void
WrapPostprocessors(py::module& post)
{
  // clang-format off
  // Volume post-processor
  auto vp = py::class_<VolumePostprocessor, std::shared_ptr<VolumePostprocessor>>(
    post,
    "VolumePostprocessor",
    R"(
    Volume post-processor.

    Wrapper of :cpp:class:`opensn::VolumePostprocessor`.
    )"
  );
  vp.def(
    py::init(
      [](py::kwargs& params)
      {
        return VolumePostprocessor::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a volume post processor object.

    Parameters
    ----------
    problem : LBSProblem
        A handle to an existing LBS problem.
    value_type : str, optional
        Type of value to compute: 'integral' (default), 'max', 'min', or 'avg'.
    )"
  );
  vp.def(
    "Execute",
    [](VolumePostprocessor& self)
    {
      self.Execute();
    },
    R"(
    Execute the postprocessor
    )"
  );
  vp.def(
    "GetValue",
    [](VolumePostprocessor& self)
    {
      const auto& arr = self.GetValue();
      auto dims = arr.dimension();

      return py::array_t<double>(
        {dims[0], dims[1]},
        {dims[1] * sizeof(double), sizeof(double)},
        arr.data(),
        py::cast(self)
      );
    },
    R"(
    Returns the value of the postprocessor.

    Rows correspond to the spatial restriction (i.e. logical volumes, if specified)
    Columns correspond to the energy restrictions (i.e. groups, or groups within a groupset if specified)
    )"
  );

  auto xs_sens =
    py::class_<CrossSectionSensitivityPostprocessor,
               std::shared_ptr<CrossSectionSensitivityPostprocessor>>(
      post,
      "CrossSectionSensitivityPostprocessor",
      R"(
    Cross-section sensitivity postprocessor.
    )"
    );
  xs_sens.def(
    py::init(
      [](py::kwargs& params)
      {
        return CrossSectionSensitivityPostprocessor::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a cross-section sensitivity postprocessor.

    sigma_t sensitivities use angular fluxes. Scattering and production
    sensitivities use flux moments.
    )"
  );
  xs_sens.def(
    "Execute",
    [](CrossSectionSensitivityPostprocessor& self)
    {
      self.Execute();
    },
    R"(
    Execute the postprocessor.
    )"
  );
  xs_sens.def(
    "GetValue",
    [](CrossSectionSensitivityPostprocessor& self)
    {
      const auto& arr = self.GetValue();
      auto dims = arr.dimension();

      return py::array_t<double>(
        {dims[0], dims[1]},
        {dims[1] * sizeof(double), sizeof(double)},
        arr.data(),
        py::cast(self)
      );
    },
    R"(
    Returns the sensitivity values.
    )"
  );
  // clang-format on
}

// Wrap the post-processing components of OpenSn.
void
py_post(py::module& pyopensn)
{
  py::module post = pyopensn.def_submodule("post", "Post-processing module.");
  WrapPostprocessors(post);
}

} // namespace opensn
