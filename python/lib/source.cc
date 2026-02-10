// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/volumetric_source/volumetric_source.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/point_source/point_source.h"
#include "framework/math/functions/function.h"
#include <memory>

namespace opensn
{

namespace
{
std::shared_ptr<GroupTimeFunction>
MakeGroupTimeFunction(const py::function& func)
{
  auto wrapper = [func](unsigned int group, double time)
  {
    py::gil_scoped_acquire gil;
    try
    {
      return func(group, time).cast<double>();
    }
    catch (const py::error_already_set& err)
    {
      if (err.matches(PyExc_TypeError))
      {
        PyErr_Clear();
        return func(group).cast<double>();
      }
      throw;
    }
  };
  return std::make_shared<GroupTimeFunction>(wrapper);
}
} // namespace

// Wrap point source
void
WrapPointSource(py::module& src)
{
  // clang-format off
  // point source
  auto point_source = py::class_<PointSource, std::shared_ptr<PointSource>>(
    src,
    "PointSource",
    R"(
    Point sources, defined by its location and a group-wise strength vector.

    Wrapper of :cpp:class:`opensn::PointSource`.
    )"
  );
  point_source.def(
    py::init(
      [](py::kwargs& params)
      {
        ParameterBlock main;
        for (auto [key, value] : params)
        {
          if (key.cast<std::string>().empty())
            continue;

          const auto name = key.cast<std::string>();
          if (name == "strength_function" and py::isinstance<py::function>(value))
          {
            main.AddParameter(name, MakeGroupTimeFunction(value.cast<py::function>()));
            continue;
          }
          main.AddParameter(pyobj_to_param_block(name, value.cast<py::object>()));
        }
        return PointSource::Create(main);
      }
    ),
    R"(
    Construct a point source from its location and strength.

    Parameters
    ----------
    location: Tuple[float, float, float]
        Coordinates of the point source.
    strength: List[float]
        Group-wise point source strength.
    strength_function: Callable
        Callable that returns the source strength for a given group (and optionally time). It must
        accept either `(group, time)` for transient problems or just `(group)` for steady-state.
        If strength_function is provided, do not specify start_time/end_time; implement any time
        dependence inside the callback instead.
    )"
  );
  // clang-format on
}

// Wrap volumetric source
void
WrapVolumetricSource(py::module& src)
{
  // clang-format off
  // volumetric source
  auto volumetric_source = py::class_<VolumetricSource, std::shared_ptr<VolumetricSource>>(
    src,
    "VolumetricSource",
    R"(
      Multi-group isotropic volumetric sources.

      Wrapper of :cpp:class:`opensn::VolumetricSource`.
    )"
  );
  volumetric_source.def(
    py::init(
      [](py::kwargs& params)
      {
        ParameterBlock main;
        for (auto [key, value] : params)
        {
          if (key.cast<std::string>().empty())
            continue;

          const auto name = key.cast<std::string>();
          if (name == "strength_function" and py::isinstance<py::function>(value))
          {
            main.AddParameter(name, MakeGroupTimeFunction(value.cast<py::function>()));
            continue;
          }
          main.AddParameter(pyobj_to_param_block(name, value.cast<py::object>()));
        }
        return VolumetricSource::Create(main);
      }
    ),
    R"(
    Construct a multi-group isotropic volumetric sources.

    Parameters
    ----------
    block_ids: List[int]
        An array of block IDs the volumetric source is present within.
    logical_volume: pyopensn.logvol.LogicalVolume
        Logical volume that the volumetric source is defined within.
    group_strength: List[float]
        An array of multi-group source strength values. Note that this is only used when a function
        is not provided.
    func: pyopensn.math.VectorSpatialFunction
        Function to be used to define the source.
    strength_function: Callable
        Callable that returns the source strength for a given group (and optionally time). It must
        accept either `(group, time)` for transient problems or just `(group)` for steady-state.
        If strength_function is provided, do not specify start_time/end_time; implement any time
        dependence inside the callback instead.
    )"
  );
  // clang-format on
}

// Wrap the source components of OpenSn
void
py_source(py::module& pyopensn)
{
  py::module src = pyopensn.def_submodule("source", "Source module.");
  WrapPointSource(src);
  WrapVolumetricSource(src);
}

} // namespace opensn
