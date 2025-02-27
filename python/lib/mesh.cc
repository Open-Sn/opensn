// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/mesh/mesh_generator/extruder_mesh_generator.h"
#include "framework/mesh/mesh_generator/orthogonal_mesh_generator.h"
#include "framework/mesh/mesh_generator/from_file_mesh_generator.h"
#include "framework/mesh/mesh_generator/split_file_mesh_generator.h"
#include "framework/mesh/mesh_generator/distributed_mesh_generator.h"
#include "framework/graphs/graph_partitioner.h"
#include "framework/graphs/kba_graph_partitioner.h"
#include "framework/graphs/linear_graph_partitioner.h"
#include "framework/graphs/petsc_graph_partitioner.h"
#include <memory>

namespace opensn
{

// clang-format off

void wrap_mesh(py::module &mesh)
{
  auto mesh_continuum = py::class_<MeshContinuum, std::shared_ptr<MeshContinuum>>(
    mesh, "MeshContinuum", R"(
    Mesh continuum.

    Wrapper of :cpp:class:`opensn::MeshContinuum`.
    )");

  mesh_continuum.def_property("dimension", &MeshContinuum::GetDimension, &MeshContinuum::SetDimension,
                              "Number of dimensions of the mesh.");
  mesh_continuum.def("SetUniformMaterialID", &MeshContinuum::SetUniformMaterialID,
                     "Set material ID's for all cells to the specified material ID.");
  mesh_continuum.def("SetMaterialIDFromLogical", &MeshContinuum::SetMaterialIDFromLogical,
                     "Set material ID's using a logical volume.");
  mesh_continuum.def("SetBoundaryIDFromLogical", &MeshContinuum::SetBoundaryIDFromLogical,
                     "Set boundary ID's using a logical volume.");
  mesh_continuum.def("SetupOrthogonalBoundaries", &MeshContinuum::SetupOrthogonalBoundaries, "...");

  auto surface_mesh = py::class_<SurfaceMesh, std::shared_ptr<SurfaceMesh>>(
    mesh, "SurfaceMesh", R"(
    Surface mesh.

    Wrapper of :cpp:class:`opensn::SurfaceMesh`.
    )");
}

void wrap_mesh_generator(py::module &mesh)
{
  auto mesh_generator = py::class_<MeshGenerator, std::shared_ptr<MeshGenerator>>(
    mesh, "MeshGenerator", R"(
    Generic mesh generator.

    Wrapper of :cpp:class:`opensn::MeshGenerator`.
    )");
  
  mesh_generator.def("Execute", &MeshGenerator::Execute, "Final execution step.");

  auto extruder_mesh_generator = py::class_<ExtruderMeshGenerator, std::shared_ptr<ExtruderMeshGenerator>, MeshGenerator>(
    mesh, "ExtruderMeshGenerator", R"(
    Extruded mesh generator.

    Wrapper of :cpp:class:`opensn::ExtruderMeshGenerator`.
    )");
  
  extruder_mesh_generator.def(py::init([](py::kwargs &params)
  {
    return ExtruderMeshGenerator::Create(kwargs_to_param_block(params));
  }), "Construct an extruded mesh generator.");

  auto orthogonal_mesh_generator = py::class_<OrthogonalMeshGenerator, std::shared_ptr<OrthogonalMeshGenerator>, MeshGenerator>(
    mesh, "OrthogonalMeshGenerator", R"(
    Orthogonal mesh generator.

    Wrapper of :cpp:class:`opensn::OrthogonalMeshGenerator`.
    )");

  orthogonal_mesh_generator.def(py::init([](py::kwargs &params)
  {
    return OrthogonalMeshGenerator::Create(kwargs_to_param_block(params));
  }), "Construct an orthogonal mesh generator.");
}

void wrap_graph_partitioner(py::module &mesh)
{
  auto graph_partitioner = py::class_<GraphPartitioner, std::shared_ptr<GraphPartitioner>>(
    mesh, "GraphPartitioner", R"(
    Generic graph partitioner.

    Wrapper of :cpp:class:`opensn::GraphPartitioner`.
    )");

  auto kba_graph_partitioner = py::class_<KBAGraphPartitioner, std::shared_ptr<KBAGraphPartitioner>, GraphPartitioner>(
    mesh, "KBAGraphPartitioner", R"(
    KBA graph partitioner.

    Wrapper of :cpp:class:`opensn::KBAGraphPartitioner`.
    )");
  
  kba_graph_partitioner.def(py::init([](py::kwargs &params)
  {
    return KBAGraphPartitioner::Create(kwargs_to_param_block(params));
  }), "Construct a KBA graph partitioner.");
}

void py_mesh(py::module &pyopensn)
{
  py::module mesh = pyopensn.def_submodule("mesh", "Mesh generation module.");
  wrap_mesh(mesh);
  wrap_mesh_generator(mesh);
  wrap_graph_partitioner(mesh);
}

// clang-format on

} // namespace opensn
