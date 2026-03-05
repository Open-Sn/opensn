// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/volume_postprocessor.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/surface_postprocessor.h"
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

  // Surface post-processor value type enum
  py::enum_<SurfacePostprocessor::ValueType>(post, "SurfacePostprocessorValueType")
    .value("INTEGRAL", SurfacePostprocessor::ValueType::INTEGRAL)
    .value("MAX", SurfacePostprocessor::ValueType::MAX)
    .value("MIN", SurfacePostprocessor::ValueType::MIN)
    .value("AVERAGE", SurfacePostprocessor::ValueType::AVERAGE);

  py::enum_<SurfacePostprocessor::CurrentType>(post, "SurfacePostprocessorCurrentType")
    .value("INCOMING", SurfacePostprocessor::CurrentType::INCOMING)
    .value("OUTGOING", SurfacePostprocessor::CurrentType::OUTGOING)
    .value("NET", SurfacePostprocessor::CurrentType::NET);

  // Surface post-processor
  auto sp = py::class_<SurfacePostprocessor, std::shared_ptr<SurfacePostprocessor>>(
    post,
    "SurfacePostprocessor",
    R"(
    Surface post-processor.

    Wrapper of :cpp:class:`opensn::SurfacePostprocessor`.
    )"
  );
  sp.def(
    py::init(
      [](py::kwargs& params)
      {
        return SurfacePostprocessor::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a surface post-processor object.

    Parameters
    ----------
    problem : DiscreteOrdinatesProblem
        A handle to an existing discrete ordinates problem.
    value_type : str, required
        Type of value to compute: 'integral', 'max', 'min', or 'avg'.
    current_type : str, required
        Type of value to compute: 'incoming', 'outgoing', or 'net'.
    )"
  );
  sp.def(
    "Execute",
    [](SurfacePostprocessor& self){
      self.Execute();
    },
    R"(
    Execute the postprocessor
    )"
  );
  sp.def(
    "GetValue",
    [](SurfacePostprocessor& self)
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
