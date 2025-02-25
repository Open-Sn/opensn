/*
 * Created on Sat, February 22
 *
 * Copyright (c) 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
 */

#include "pyapi.hpp"

#include <memory>     // std::shared_ptr
#include <stdexcept>  // std::invalid_argument
#include <string>     // std::string, std::to_string

#include "framework/graphs/graph_partitioner.h"                                               // opensn::GraphPartitioner
#include "framework/materials/material_property.h"                                            // opensn::MaterialProperty
#include "framework/math/quadratures/angular/angular_quadrature.h"                            // opensn::AngularQuadrature
#include "framework/mesh/logical_volume/logical_volume.h"                                     // opensn::LogicalVolume
#include "framework/mesh/mesh_continuum/mesh_continuum.h"                                     // opensn::MeshContinuum
#include "framework/mesh/mesh_generator/mesh_generator.h"                                     // opensn::MeshGenerator
#include "framework/mesh/surface_mesh/surface_mesh.h"                                         // opensn::SurfaceMesh
#include "framework/physics/solver.h"                                                         // opensn::Solver
#include "modules/linear_boltzmann_solvers/lbs_solver/point_source/point_source.h"            // opensn::PointSource
#include "modules/linear_boltzmann_solvers/lbs_solver/volumetric_source/volumetric_source.h"  // opensn::VolumetricSource

#define TO_PARAMBLOCK(class_name)                                                                                      \
    if (py::isinstance<class_name>(obj)) return ParameterBlock(key, obj.cast<std::shared_ptr<class_name>>())           \

namespace opensn {

// Convert a Python obejct into a ParameterBlock
ParameterBlock pyobj_to_param_block(const std::string & key, const py::object & obj) {
    // basic types
    if (py::isinstance<py::bool_>(obj)) {
        return ParameterBlock(key, obj.cast<bool>());
    }
    if (py::isinstance<py::int_>(obj)) {
        return ParameterBlock(key, obj.cast<int>());
    }
    if (py::isinstance<py::float_>(obj)) {
        return ParameterBlock(key, obj.cast<double>());
    }
    if (py::isinstance<py::str>(obj)) {
        return ParameterBlock(key, obj.cast<std::string>());
    }
    // dictionary (recursive add)
    if (py::isinstance<py::dict>(obj)) {
        ParameterBlock main(key);
        py::dict dict_obj = obj.cast<py::dict>();
        for (auto [key, value] : dict_obj) {
            main.AddParameter(pyobj_to_param_block(key.cast<std::string>(), value.cast<py::object>()));
        }
        return main;
    }
    // list or tuple (recursive add)
    if (py::isinstance<py::list>(obj) || py::isinstance<py::tuple>(obj)) {
        ParameterBlock list(key);
        list.ChangeToArray();
        for (py::handle element : obj) {
            std::string index_string = std::to_string(list.GetNumParameters());
            list.AddParameter(pyobj_to_param_block(index_string, element.cast<py::object>()));
        }
        return list;
    }
    // check for each class in the OpenSn library
    TO_PARAMBLOCK(AngularQuadrature);
    TO_PARAMBLOCK(GraphPartitioner);
    TO_PARAMBLOCK(LogicalVolume);
    TO_PARAMBLOCK(MaterialProperty);
    TO_PARAMBLOCK(MeshContinuum);
    TO_PARAMBLOCK(MeshGenerator);
    TO_PARAMBLOCK(SurfaceMesh);
    TO_PARAMBLOCK(PointSource);
    TO_PARAMBLOCK(VolumetricSource);
    TO_PARAMBLOCK(Solver);
    // throw and return
    throw std::invalid_argument("Unsupported argument type.");
    return ParameterBlock();
}

// Translate a Python dictionary into a ParameterBlock
ParameterBlock kwargs_to_param_block(const py::kwargs & params) {
    // initialize main parameter dict
    ParameterBlock main;
    // append corresponding parameters to the main
    for (auto [key, value] : params) {
        // skip empty keys
        if (key.cast<std::string>().empty()) {
            continue;
        }
        // recursively add parameters
        main.AddParameter(pyobj_to_param_block(key.cast<std::string>(), value.cast<py::object>()));
    }
    return main;
}

}  // namespace opensn
