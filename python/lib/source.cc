// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/volumetric_source/volumetric_source.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/point_source/point_source.h"
#include <memory>

namespace opensn
{

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
        return PointSource::Create(kwargs_to_param_block(params));
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
        return VolumetricSource::Create(kwargs_to_param_block(params));
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
    group_strength: List[double]
        An array of multi-group source strength values. Note that this is only used when a function
        is not provided.
    func: pyopensn.math.VectorSpatialFunction
        Function to be used to define the source.
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
