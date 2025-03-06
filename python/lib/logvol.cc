// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/logical_volume/boolean_logical_volume.h"
#include "framework/mesh/logical_volume/rcc_logical_volume.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/mesh/logical_volume/sphere_logical_volume.h"
#include "framework/mesh/logical_volume/surface_mesh_logical_volume.h"
#include <memory>

namespace opensn
{

// Wrap logical volume
void WrapLogicalVolume(py::module& logvol)
{
  // generic logical volume
  auto logical_volume = py::class_<LogicalVolume, std::shared_ptr<LogicalVolume>>(
    logvol,
    "LogicalVolume",
    R"(
    Generic logical volume.

    Wrapper of :cpp:class:`opensn::LogicalVolume`.
    )"
  );
  logical_volume.def(
    "Inside",
    &LogicalVolume::Inside,
    "Check if a point is inside or outside the logical volume.",
    py::arg("point")
  );

  // boolean logical volume
  auto boolean_logical_volume = py::class_<BooleanLogicalVolume,
                                           std::shared_ptr<BooleanLogicalVolume>, LogicalVolume>(
    logvol,
    "BooleanLogicalVolume",
    R"(
    Boolean logical volume.

    Wrapper of :cpp:class:`opensn::BooleanLogicalVolume`.
    )"
  );
  boolean_logical_volume.def(
    py::init(
      [](py::kwargs& params)
      {
        return BooleanLogicalVolume::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a boolean logical volume from a list of combinatorial logics. Each combinatorial logic
    are represented by a Python dictionary with two keys:

     - ``op`` maps to a boolean. True means included and False means excluded.
     - ``lv`` maps to a :py:class:`pyopensn.logvol.LogicalVolume`.

    Parameters
    ----------
    parts: List[Dict]
        List of dictionary of combinatorial logics.

    Examples
    --------
    >>> blv = BooleanLogicalVolume(
    ...     parts=[{ "op": True, "lv": volume1 }, { "op": False, "lv": volume2 }]
    ... )
    )"
  );

  // RCC logical volume
  auto rcc_logical_volume = py::class_<RCCLogicalVolume, std::shared_ptr<RCCLogicalVolume>,
                                       LogicalVolume>(
    logvol,
    "RCCLogicalVolume",
    R"(
    Right circular cylinder logical volume.

    Wrapper of :cpp:class:`opensn::RCCLogicalVolume`.
    )"
  );
  rcc_logical_volume.def(
    py::init(
      [](py::kwargs& params)
      {
        return RCCLogicalVolume::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct an RCC logical volume.

    Parameters
    ----------
    r: float, default=1.0
        Radius of the sphere.
    x0: float, default=0.0
        X-coordinate of the volume base.
    y0: float, default=0.0
        Y-coordinate of the volume base.
    z0: float, default=0.0
        Z-coordinate of the volume base.
    vx: float, default=0.0
        X-component of the volume extrusion vector.
    vy: float, default=0.0
        Y-component of the volume extrusion vector.
    vz: float, default=1.0
        Z-component of the volume extrusion vector.
    )"
  );

  // RPP logical volume
  auto rpp_logical_volume = py::class_<RPPLogicalVolume, std::shared_ptr<RPPLogicalVolume>,
                                       LogicalVolume>(
    logvol,
    "RPPLogicalVolume",
    R"(
    Rectangular parallel piped logical volume.

    Wrapper of :cpp:class:`opensn::RPPLogicalVolume`.
    )"
  );
  rpp_logical_volume.def(
    py::init(
      [](py::kwargs& params)
      {
        return RPPLogicalVolume::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct an RPP logical volume.

    Parameters
    ----------
    xmin: float, default=0.0
        X-min of the volume.
    xmax: float, default=1.0
        X-max of the volume.
    ymin: float, default=0.0
        Y-min of the volume.
    ymax: float, default=1.0
        Y-max of the volume.
    zmin: float, default=0.0
        Z-min of the volume.
    zmax: float, default=1.0
        Z-max of the volume.
    infx: bool, default=False
        Flag, when true, will ignore xmin and xmax.
    infy: bool, default=False
        Flag, when true, will ignore ymin and ymax.
    infz: bool, default=False
        Flag, when true, will ignore zmin and zmax.
    )"
  );

  // spherical logical volume
  auto spherical_logical_volume = py::class_<SphereLogicalVolume,
                                             std::shared_ptr<SphereLogicalVolume>, LogicalVolume>(
    logvol,
    "SphereLogicalVolume",
    R"(
    Spherical logical volume.

    Wrapper of :cpp:class:`opensn::SphereLogicalVolume`.
    )"
  );
  spherical_logical_volume.def(
    py::init(
      [](py::kwargs& params)
      {
        return SphereLogicalVolume::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a spherical logical volume.

    Parameters
    ----------
    r: float, default=1.0
        Radius of the sphere.
    x: float, default=0.0
        X-location of the volume.
    y: float, default=0.0
        Y-location of the volume.
    z: float, default=0.0
        Z-location of the volume.
    )"
  );

  // surface mesh logical volume
  auto surface_mesh_logical_volume = py::class_<SurfaceMeshLogicalVolume,
                                                std::shared_ptr<SurfaceMeshLogicalVolume>,
                                                LogicalVolume>(
    logvol,
    "SurfaceMeshLogicalVolume",
    R"(
    Surface mesh logical volume.

    Wrapper of :cpp:class:`opensn::SurfaceMeshLogicalVolume`.
    )"
  );
  surface_mesh_logical_volume.def(
    py::init(
      [](py::kwargs& params)
      {
        return SurfaceMeshLogicalVolume::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a surface mesh logical volume.

    Parameters
    ----------
    surface_mesh: pyopensn.mesh.SurfaceMesh
        Associated surface mesh.
    )"
  );
}

// Wrap the logical volume components of OpenSn
void py_logvol(py::module &pyopensn)
{
  py::module logvol = pyopensn.def_submodule("logvol", "Logical volume module.");
  WrapLogicalVolume(logvol);
}

} // namespace opensn
