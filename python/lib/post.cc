// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/volume_postprocessor.h"
#include <pybind11/stl.h>

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
      return self.GetValue();
    },
    R"(
    Returns the value of the postprocessor.

    Rows correspond to the spatial restriction (i.e. logical volumes, if specified)
    Columns correspond to the energy restrictions (i.e. groups, or groups within a groupset if specified)
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
