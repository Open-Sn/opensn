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
  // Volume post-processor value type enum
  py::enum_<VolumePostprocessor::ValueType>(post, "VolumePostprocessorValueType")
    .value("INTEGRAL", VolumePostprocessor::ValueType::INTEGRAL)
    .value("MAX", VolumePostprocessor::ValueType::MAX)
    .value("MIN", VolumePostprocessor::ValueType::MIN)
    .value("AVERAGE", VolumePostprocessor::ValueType::AVERAGE);

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
    "Initialize",
    [](VolumePostprocessor& self){
      self.Initialize();
    },
    R"(
      TODO: finish this
    )"
  );
  vp.def(
    "Execute",
    [](VolumePostprocessor& self){
      self.Execute();
    },
    R"(
      TODO: finish this
    )"
  );
  vp.def(
    "GetValue",
    [](VolumePostprocessor& self)
    {
      return self.GetValue();
    },
    R"(
    TODO: finish this
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
