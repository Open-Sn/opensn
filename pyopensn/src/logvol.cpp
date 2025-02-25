// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "pyapi.hpp"

#include <memory>

#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/logical_volume/boolean_logical_volume.h"
#include "framework/mesh/logical_volume/rcc_logical_volume.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/mesh/logical_volume/sphere_logical_volume.h"
#include "framework/mesh/logical_volume/surface_mesh_logical_volume.h"

namespace opensn {

// Wrap logical volume
static void wrap_logical_volume(py::module & logvol) {
    // generic logical volume
    auto logical_volume = py::class_<LogicalVolume, std::shared_ptr<LogicalVolume>>(
        logvol,
        "LogicalVolume",
        R"(
        Generic logical volume.

        Wrapper of :cpp:class:`opensn::LogicalVolume`.
        )"
    );
    // boolean logical volume
    auto boolean_logical_volume = py::class_<BooleanLogicalVolume, std::shared_ptr<BooleanLogicalVolume>, LogicalVolume>(
        logvol,
        "BooleanLogicalVolume",
        R"(
        Boolean logical volume.

        Wrapper of :cpp:class:`opensn::BooleanLogicalVolume`.
        )"
    );
    // RCC logical volume
    auto rcc_logical_volume = py::class_<RCCLogicalVolume, std::shared_ptr<RCCLogicalVolume>, LogicalVolume>(
        logvol,
        "RCCLogicalVolume",
        R"(
        Right circular cylinder logical volume.

        Wrapper of :cpp:class:`opensn::RCCLogicalVolume`.
        )"
    );
    // RPP logical volume
    auto rpp_logical_volume = py::class_<RPPLogicalVolume, std::shared_ptr<RPPLogicalVolume>, LogicalVolume>(
        logvol,
        "RPPLogicalVolume",
        R"(
        Rectangular parallel piped logical volume.

        Wrapper of :cpp:class:`opensn::RPPLogicalVolume`.
        )"
    );
    // Spherical logical volume
    auto spherical_logical_volume = py::class_<SphereLogicalVolume, std::shared_ptr<SphereLogicalVolume>, LogicalVolume>(
        logvol,
        "SphereLogicalVolume",
        R"(
        Spherical logical volume.

        Wrapper of :cpp:class:`opensn::SphereLogicalVolume`.
        )"
    );
    // Surface mesh logical volume
    auto surface_mesh_logical_volume = py::class_<SurfaceMeshLogicalVolume, std::shared_ptr<SurfaceMeshLogicalVolume>, LogicalVolume>(
        logvol,
        "SurfaceMeshLogicalVolume",
        R"(
        Surface mesh logical volume.

        Wrapper of :cpp:class:`opensn::SurfaceMeshLogicalVolume`.
        )"
    );
}

// Wrap the logical volume components of OpenSn
void py_logvol(py::module & pyopensn) {
    py::module logvol = pyopensn.def_submodule("logvol", "Logical volume module.");
    wrap_logical_volume(logvol);
}

}  // namespace opensn
