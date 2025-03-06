// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"

#include <memory>

#include "framework/graphs/graph_partitioner.h"
#include "framework/graphs/kba_graph_partitioner.h"
#include "framework/graphs/linear_graph_partitioner.h"
#include "framework/graphs/petsc_graph_partitioner.h"
#include "framework/mesh/io/mesh_io.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/mesh/mesh_generator/extruder_mesh_generator.h"
#include "framework/mesh/mesh_generator/orthogonal_mesh_generator.h"
#include "framework/mesh/mesh_generator/from_file_mesh_generator.h"
#include "framework/mesh/mesh_generator/split_file_mesh_generator.h"
#include "framework/mesh/mesh_generator/distributed_mesh_generator.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"

namespace opensn
{

// Wrap mesh continuum
void WrapMesh(py::module& mesh)
{
  // Mesh continuum
  auto mesh_continuum = py::class_<MeshContinuum, std::shared_ptr<MeshContinuum>>(
    mesh,
    "MeshContinuum",
    R"(
    Mesh continuum.

    Wrapper of :cpp:class:`opensn::MeshContinuum`.
    )"
  );
  mesh_continuum.def_property(
    "dimension",
    &MeshContinuum::GetDimension,
    &MeshContinuum::SetDimension,
    "Number of dimensions of the mesh."
  );
  mesh_continuum.def(
    "SetUniformMaterialID",
    &MeshContinuum::SetUniformMaterialID,
    "Set material ID's for all cells to the specified material ID.",
    py::arg("mat_id")
  );
  mesh_continuum.def(
    "SetMaterialIDFromLogicalVolume",
    &MeshContinuum::SetMaterialIDFromLogical,
    R"(
    Set material ID's using a logical volume.

    Parameters
    ----------
    log_vol: pyopensn.logvol.LogicalVolume
        ???
    mat_id: int
        ???
    sense: bool
        ???
    )",
    py::arg("log_vol"),
    py::arg("mat_id"),
    py::arg("sense")
  );
  mesh_continuum.def(
    "SetBoundaryIDFromLogicalVolume",
    &MeshContinuum::SetBoundaryIDFromLogical,
    R"(
    Set boundary ID's using a logical volume.

    Parameters
    ----------

    log_vol: pyopensn.logvol.LogicalVolume
        ???
    boundary_name: str
        ???
    sense: bool
        ???
    )",
    py::arg("log_vol"),
    py::arg("boundary_name"),
    py::arg("sense") = true
  );
  mesh_continuum.def(
    "SetupOrthogonalBoundaries",
    &MeshContinuum::SetupOrthogonalBoundaries,
    "???"
  );
  mesh_continuum.def(
    "ExportToPVTU",
    [](std::shared_ptr<MeshContinuum> self, const std::string& file_name) {
        MeshIO::ToPVTU(self, file_name);
    },
    R"(
    Write grid cells into PVTU format.

    Parameters
    ----------
    file_base_name: str
        Base name of the output file.
    )",
    py::arg("file_base_name")
  );

  // Surface mesh
  auto surface_mesh = py::class_<SurfaceMesh, std::shared_ptr<SurfaceMesh>>(
    mesh,
    "SurfaceMesh",
    R"(
    Surface mesh.

    Wrapper of :cpp:class:`opensn::SurfaceMesh`.
    )"
  );
  surface_mesh.def(
    py::init(
      [](py::kwargs& params)
      {
        return SurfaceMesh::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a surface mesh.

    This function does not take any arguments.
    )"
  );
  surface_mesh.def(
    "ImportFromOBJFile",
    &SurfaceMesh::ImportFromOBJFile,
    R"(
    Loads a surface mesh from a wavefront .obj file.

    Parameters
    ----------
    file_name: str
        ???
    as_poly: bool, default=False
        ???
    transform: pyopensn.math.Vector3, default=(0.0, 0.0, 0.0)
        ???
    )",
    py::arg("file_name"),
    py::arg("as_poly") = false,
    py::arg("transform") = Vector3()
  );
  surface_mesh.def(
    "ImportFromTriangleFiles",
    &SurfaceMesh::ImportFromTriangleFiles,
    R"(
    Loads a surface mesh from triangle's file format.

    Parameters
    ----------
    file_name: str
        ???
    as_poly: bool
        ???
    )",
    py::arg("file_name"),
    py::arg("as_poly")
  );
  surface_mesh.def(
    "ImportFromMshFiles",
    &SurfaceMesh::ImportFromMshFiles,
    R"(
    Loads a surface mesh from gmsh's file format.

    Parameters
    ----------
    file_name: str
        ???
    as_poly: bool
        ???
    )",
    py::arg("file_name"),
    py::arg("as_poly")
  );
}

// Wrap mesh generator
void WrapMeshGenerator(py::module& mesh)
{
  // Base mesh generator
  auto mesh_generator = py::class_<MeshGenerator, std::shared_ptr<MeshGenerator>>(
    mesh,
    "MeshGenerator",
    R"(
    Generic mesh generator.

    Wrapper of :cpp:class:`opensn::MeshGenerator`.
    )"
  );
  mesh_generator.def(
    "Execute",
    &MeshGenerator::Execute,
    "Final execution step."
  );

  // Extruded mesh generator
  auto extruder_mesh_generator = py::class_<ExtruderMeshGenerator,
                                            std::shared_ptr<ExtruderMeshGenerator>, MeshGenerator>(
    mesh,
    "ExtruderMeshGenerator",
    R"(
    Extruded mesh generator.

    Wrapper of :cpp:class:`opensn::ExtruderMeshGenerator`.
    )"
  );
  extruder_mesh_generator.def(
    py::init(
      [](py::kwargs& params)
      {
        return ExtruderMeshGenerator::Create(kwargs_to_param_block(params));
      }
    ),
    "Construct an extruded mesh generator."
  );

  // Orthogonal mesh generator
  auto orthogonal_mesh_generator = py::class_<OrthogonalMeshGenerator,
                                              std::shared_ptr<OrthogonalMeshGenerator>,
                                              MeshGenerator>(
    mesh,
    "OrthogonalMeshGenerator",
    R"(
    Orthogonal mesh generator.

    Wrapper of :cpp:class:`opensn::OrthogonalMeshGenerator`.
    )"
  );
  orthogonal_mesh_generator.def(
    py::init(
      [](py::kwargs& params)
      {
        return OrthogonalMeshGenerator::Create(kwargs_to_param_block(params));
      }
    ),
    "Construct an orthogonal mesh generator."
  );
}

// Wrap graph partitioner
void WrapGraphPartitioner(py::module& mesh)
{
  // Base graph partitioner
  auto graph_partitioner = py::class_<GraphPartitioner, std::shared_ptr<GraphPartitioner>>(
    mesh,
    "GraphPartitioner",
    R"(
    Generic graph partitioner.

    Wrapper of :cpp:class:`opensn::GraphPartitioner`.
    )");

  // KBA graph partitioner
  auto kba_graph_partitioner = py::class_<KBAGraphPartitioner, std::shared_ptr<KBAGraphPartitioner>,
                                          GraphPartitioner>(
    mesh,
    "KBAGraphPartitioner",
    R"(
    KBA graph partitioner.

    Wrapper of :cpp:class:`opensn::KBAGraphPartitioner`.
    )"
  );
  kba_graph_partitioner.def(
    py::init(
      [](py::kwargs& params)
      {
        return KBAGraphPartitioner::Create(kwargs_to_param_block(params));
      }
    ),
    "Construct a KBA graph partitioner."
  );

  // PETSc graph partitioner
  auto petsc_graph_partitioner = py::class_<PETScGraphPartitioner,
                                            std::shared_ptr<PETScGraphPartitioner>,
                                            GraphPartitioner>(
    mesh,
    "PETScGraphPartitioner",
    R"(
    PETSc graph partitioner.

    Wrapper of :cpp:class:`opensn::PETScGraphPartitioner`.
    )"
  );
  petsc_graph_partitioner.def(
    py::init(
      [](py::kwargs & params) {
        return PETScGraphPartitioner::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a linear graph partitioner.
    )"
  );
}

void py_mesh(py::module& pyopensn)
{
  py::module mesh = pyopensn.def_submodule("mesh", "Mesh generation module.");
  WrapMesh(mesh);
  WrapMeshGenerator(mesh);
  WrapGraphPartitioner(mesh);
}

} // namespace opensn
