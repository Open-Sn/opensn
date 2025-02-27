// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/volumetric_source/volumetric_source.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/point_source/point_source.h"
#include <memory>

namespace opensn
{

// clang-format off

// Wrap source
void wrap_source(py::module &src)
{
  // point source
  auto point_source = py::class_<PointSource, std::shared_ptr<PointSource>>(
    src, "PointSource",
    R"(
    Point sources, defined by its location and a group-wise strength vector.

    Wrapper of :cpp:class:`opensn::PointSource`.
    )");

  point_source.def(
    py::init(
      [](py::kwargs &params)
      {
        return PointSource::Create(kwargs_to_param_block(params));
      }),
    R"(
    Construct a point source from its location and strength.

    Parameters
    ----------
    location: Tuple[float, float, float]
        Coordinates of the point source.
    strength: List[float]
        Group-wise point source strength.
    )");

  // volumetric source
  auto volumetric_source = py::class_<VolumetricSource, std::shared_ptr<VolumetricSource>>(
    src, "VolumetricSource",
    R"(
    Multi-group isotropic volumetric sources.

    Wrapper of :cpp:class:`opensn::VolumetricSource`.
    )");

  volumetric_source.def(
    py::init(
      [](py::kwargs &params)
      {
        return VolumetricSource::Create(kwargs_to_param_block(params));
      }),
    R"(
    Construct a multi-group isotropic volumetric sources.

    Parameters
    ----------
    block_ids: List[int]
        ???
    logical_volume: opensn.LogicalVolume
        ???
    group_strength: List[double]
        ???
    func: ???
        ???
    )");
}

// Wrap the source components of OpenSn
void py_source(py::module &pyopensn)
{
  py::module src = pyopensn.def_submodule("source", "Source module.");
  wrap_source(src);
}

// clang-format on

} // namespace opensn
