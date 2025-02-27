// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "py_api.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/point_source/point_source.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/volumetric_source/volumetric_source.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/graphs/graph_partitioner.h"
#include "framework/materials/material_property.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/physics/solver.h"
#include <memory>
#include <stdexcept>
#include <string>

#define TO_PARAMBLOCK(class_name)                                                                  \
  if (py::isinstance<class_name>(obj))                                                             \
  return ParameterBlock(key, obj.cast<std::shared_ptr<class_name>>())

namespace opensn
{

// Convert a Python object into a ParameterBlock
ParameterBlock
pyobj_to_param_block(const std::string& key, const py::object& obj)
{
  // basic types
  if (py::isinstance<py::bool_>(obj))
  {
    return ParameterBlock(key, obj.cast<bool>());
  }
  if (py::isinstance<py::int_>(obj))
  {
    return ParameterBlock(key, obj.cast<int>());
  }
  if (py::isinstance<py::float_>(obj))
  {
    return ParameterBlock(key, obj.cast<double>());
  }
  if (py::isinstance<py::str>(obj))
  {
    return ParameterBlock(key, obj.cast<std::string>());
  }

  // dictionary (recursive add)
  if (py::isinstance<py::dict>(obj))
  {
    ParameterBlock main(key);
    py::dict dict_obj = obj.cast<py::dict>();
    for (auto [key, value] : dict_obj)
    {
      main.AddParameter(pyobj_to_param_block(key.cast<std::string>(), value.cast<py::object>()));
    }
    return main;
  }

  // list or tuple (recursive add)
  if (py::isinstance<py::list>(obj) || py::isinstance<py::tuple>(obj))
  {
    ParameterBlock list(key);
    list.ChangeToArray();
    for (py::handle element : obj)
    {
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
ParameterBlock
kwargs_to_param_block(const py::kwargs& params)
{
  // initialize main parameter dict
  ParameterBlock main;

  // append corresponding parameters to the main
  for (auto [key, value] : params)
  {
    // skip empty keys
    if (key.cast<std::string>().empty())
    {
      continue;
    }

    // recursively add parameters
    main.AddParameter(pyobj_to_param_block(key.cast<std::string>(), value.cast<py::object>()));
  }

  return main;
}

} // namespace opensn
