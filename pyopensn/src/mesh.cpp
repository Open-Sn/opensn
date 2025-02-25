// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "pyapi.hpp"

#include <memory>

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

namespace opensn {

// Wrap mesh
static void wrap_mesh(py::module & mesh) {
    // mesh continuum
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
        "Set material ID's for all cells to the specified material ID."
    );
    mesh_continuum.def(
        "SetMaterialIDFromLogical",
        &MeshContinuum::SetMaterialIDFromLogical,
        "Set material ID's using a logical volume."
    );
    mesh_continuum.def(
        "SetBoundaryIDFromLogical",
        &MeshContinuum::SetBoundaryIDFromLogical,
        "Set boundary ID's using a logical volume."
    );
    mesh_continuum.def(
        "SetupOrthogonalBoundaries",
        &MeshContinuum::SetupOrthogonalBoundaries,
        "..."
    );
    // surface mesh
    auto surface_mesh = py::class_<SurfaceMesh, std::shared_ptr<SurfaceMesh>>(
        mesh,
        "SurfaceMesh",
        R"(
        Surface mesh.

        Wrapper of :cpp:class:`opensn::SurfaceMesh`.
        )"
    );
}

// Wrap mesh generators
static void wrap_mesh_generator(py::module & mesh) {
    // generic mesh generator
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
        R"(
        Final execution step.
        )"
    );
    // extruded mesh generator
    auto extruder_mesh_generator = py::class_<ExtruderMeshGenerator, std::shared_ptr<ExtruderMeshGenerator>, MeshGenerator>(
        mesh,
        "ExtruderMeshGenerator",
        R"(
        Extruded mesh generator.

        Wrapper of :cpp:class:`opensn::ExtruderMeshGenerator`.
        )"
    );
    extruder_mesh_generator.def(
        py::init(
            [](py::kwargs & params) {
                return ExtruderMeshGenerator::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct an extruded mesh generator.
        )"
    );
    // orthogonal mesh generator
    auto orthogonal_mesh_generator = py::class_<OrthogonalMeshGenerator, std::shared_ptr<OrthogonalMeshGenerator>, MeshGenerator>(
        mesh,
        "OrthogonalMeshGenerator",
        R"(
        Orthogonal mesh generator.

        Wrapper of :cpp:class:`opensn::OrthogonalMeshGenerator`.
        )"
    );
    orthogonal_mesh_generator.def(
        py::init(
            [](py::kwargs & params) {
                return OrthogonalMeshGenerator::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct an orthogonal mesh generator.

        Parameters
        ----------
        node_sets: List[List[float]]
            Coordinates of nodes on each axis.
        partitioner: pyopensn.GraphPartitioner
            Graph partitioning algorithm.
        )"
    );
    // from file mesh generator
    auto from_file_mesh_generator = py::class_<FromFileMeshGenerator, std::shared_ptr<FromFileMeshGenerator>, MeshGenerator>(
        mesh,
        "FromFileMeshGenerator",
        R"(
        From file mesh generator.

        Wrapper of :cpp:class:`opensn::FromFileMeshGenerator`.
        )"
    );
    from_file_mesh_generator.def(
        py::init(
            [](py::kwargs & params) {
                return FromFileMeshGenerator::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct a from-file mesh generator.
        )"
    );
    // split file mesh generator
    auto split_file_mesh_generator = py::class_<SplitFileMeshGenerator, std::shared_ptr<SplitFileMeshGenerator>, MeshGenerator>(
        mesh,
        "SplitFileMeshGenerator",
        R"(
        Split file mesh generator.

        Wrapper of :cpp:class:`opensn::SplitFileMeshGenerator`.
        )"
    );
    split_file_mesh_generator.def(
        py::init(
            [](py::kwargs & params) {
                return SplitFileMeshGenerator::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct a split-file mesh generator.
        )"
    );
    // distributed mesh generator
    auto distributed_mesh_generator = py::class_<DistributedMeshGenerator, std::shared_ptr<DistributedMeshGenerator>, MeshGenerator>(
        mesh,
        "DistributedMeshGenerator",
        R"(
        This class is responsible for generating a mesh, partitioning it, and distributing the individual partitions to
        different MPI locations. The mesh is generated on location 0, partitioned into multiple parts, serialized, and
        distributed to all other MPI ranks.

        Wrapper of :cpp:class:`opensn::DistributedMeshGenerator`.
        )"
    );
    distributed_mesh_generator.def(
        py::init(
            [](py::kwargs & params) {
                return DistributedMeshGenerator::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct a distributed mesh generator.
        )"
    );
}

// Wrap graph partitioner
static void wrap_graph_partitioner(py::module & mesh) {
    // generic graph partitioner
    auto graph_partitioner = py::class_<GraphPartitioner, std::shared_ptr<GraphPartitioner>>(
        mesh,
        "GraphPartitioner",
        R"(
        Generic graph partitioner.

        Wrapper of :cpp:class:`opensn::GraphPartitioner`.
        )"
    );
    // KBA graph partitioner
    auto kba_graph_partitioner = py::class_<KBAGraphPartitioner, std::shared_ptr<KBAGraphPartitioner>, GraphPartitioner>(
        mesh,
        "KBAGraphPartitioner",
        R"(
        KBA graph partitioner.

        Wrapper of :cpp:class:`opensn::KBAGraphPartitioner`.
        )"
    );
    kba_graph_partitioner.def(
        py::init(
            [](py::kwargs & params) {
                return KBAGraphPartitioner::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct a KBA graph partitioner.

        Parameters
        ----------
        nx: int
            Number of partitions on the x-axis.
        ny: int
            Number of partitions on the y-axis.
        xcuts: List[float]
            X-coordinates of the partition boundaries. The size of the list must be equal to (`nx' - 1).
        ycuts: List[float]
            Y-coordinates of the partition boundaries. The size of the list must be equal to (`ny' - 1).
        )"
    );
    // linear graph partitioner
    auto linear_graph_partitioner = py::class_<LinearGraphPartitioner, std::shared_ptr<LinearGraphPartitioner>, GraphPartitioner>(
        mesh,
        "LinearGraphPartitioner",
        R"(
        Linear graph partitioner.

        Wrapper of :cpp:class:`opensn::LinearGraphPartitioner`.
        )"
    );
    linear_graph_partitioner.def(
        py::init(
            [](py::kwargs & params) {
                return LinearGraphPartitioner::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct a linear graph partitioner.
        )"
    );
    // PETSc graph partitioner
    auto petsc_graph_partitioner = py::class_<PETScGraphPartitioner, std::shared_ptr<PETScGraphPartitioner>, GraphPartitioner>(
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

// Wrap the mesh components of OpenSn
void py_mesh(py::module & pyopensn) {
    py::module mesh = pyopensn.def_submodule("mesh", "Mesh generation module.");
    wrap_mesh(mesh);
    wrap_mesh_generator(mesh);
    wrap_graph_partitioner(mesh);
}

}  // namespace opensn
