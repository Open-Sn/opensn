// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "pyapi.hpp"

#include <memory>

#include "framework/materials/material_property.h"
#include "framework/materials/material.h"
#include "framework/logging/log.h"

namespace opensn {

// Wrap material
static void wrap_material(py::module & mat) {
    // material property
    auto material_property = py::class_<MaterialProperty, std::shared_ptr<MaterialProperty>>(
        mat,
        "MaterialProperty",
        R"(
        Base class for material properties.

        Wrapper of :cpp:class:`opensn::MaterialProperty`.
        )"
    );
    // material
    auto material = py::class_<Material, std::shared_ptr<Material>>(
        mat,
        "Material",
        R"(
        Base class for materials used in physics simulations.

        Wrapper of :cpp:class:`opensn::Material`.
        )"
    );
    material.def(
        "SetTransportXSections",
        &Material::SetTransportXSections,
        "Set multi-group cross sections to the material."
    );
    // add material
    mat.def(
        "AddMaterial",
        [](const std::string & name) {
            std::shared_ptr<Material> new_material = std::make_shared<Material>();
            new_material->name = name;
            material_stack.push_back(new_material);
            size_t index = material_stack.size() - 1;
            log.Log0Verbose1() << "New material added at index " << index << " with name \"" << new_material->name
                               << "\"";
            return new_material;
        },
        "Add new material with a given name."
    );
}

// Wrap the material components of OpenSn
void py_mat(py::module & pyopensn) {
    py::module mat = pyopensn.def_submodule("mat", "Material module.");
    wrap_material(mat);
}

}  // namespace opensn
