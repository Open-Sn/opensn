// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
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
#include "framework/utils/timer.h"
#include <pybind11/functional.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace opensn
{

// Wrap mesh continuum
void
WrapMesh(py::module& mesh)
{
  // clang-format off
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
    "SetUniformBlockID",
    &MeshContinuum::SetUniformBlockID,
    "Set block ID's for all cells to the specified block ID.",
    py::arg("mat_id")
  );
  mesh_continuum.def(
    "SetBlockIDFromLogical",
    &MeshContinuum::SetBlockIDFromLogical,
    R"(
    Set block ID's using a logical volume.

    Parameters
    ----------
    log_vol: pyopensn.logvol.LogicalVolume
        Logical volume that determines which mesh cells will be selected.
    mat_id: int
        Block ID that will be assigned.
    inside: bool
        If true, the selected mesh cells are the ones whose centroids are inside the logival volume.
        Otherwise, the selected meshes are the ones whose centroids are outside of the logical
        volume.
    )",
    py::arg("log_vol"),
    py::arg("mat_id"),
    py::arg("inside")
  );
  mesh_continuum.def(
    "SetBoundaryIDFromLogical",
    &MeshContinuum::SetBoundaryIDFromLogical,
    R"(
    Set boundary ID's using a logical volume.

    Parameters
    ----------
    log_vol: pyopensn.logvol.LogicalVolume
        Logical volume that determines which mesh cells will be selected.
    boundary_name: str
        Name of the boundary.
    inside: bool
        If true, the selected cell facess are the ones whose centroids are inside the logival
        volume. Otherwise, the selected meshes are the ones whose centroids are outside of the
        logical volume.
    )",
    py::arg("log_vol"),
    py::arg("boundary_name"),
    py::arg("inside") = true
  );
  mesh_continuum.def(
    "SetupOrthogonalBoundaries",
    &MeshContinuum::SetupOrthogonalBoundaries,
    "Setup boundary IDs for xmin/xmax, ymin/ymax, zmin/zmax for a right parallelpiped domain."
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
  mesh_continuum.def(
    "ComputeVolumePerBlockID",
    &MeshContinuum::ComputeVolumePerBlockID,
    "Compute volume per block ID."
  );
  mesh_continuum.def(
    "SetBlockIDFromFunction",
    [](MeshContinuum& self, const std::function<int(Vector3, int)>& func)
    {
      int local_num_cells_modified = 0;
      // change local cells
      for (Cell& cell : self.local_cells)
      {
        int new_matid = func(cell.centroid, cell.block_id);
        if (cell.block_id != new_matid)
        {
          cell.block_id = new_matid;
          ++local_num_cells_modified;
        }
      }
      // change ghost cells
      std::vector<std::uint64_t> ghost_ids = self.cells.GetGhostGlobalIDs();
      for (std::uint64_t ghost_id : ghost_ids)
      {
        Cell& cell = self.cells[ghost_id];
        int new_matid = func(cell.centroid, cell.block_id);
        if (cell.block_id != new_matid)
        {
          cell.block_id = new_matid;
          ++local_num_cells_modified;
        }
      }
      // print number of modified cells
      int global_num_cells_modified;
      mpi_comm.all_reduce(local_num_cells_modified, global_num_cells_modified, mpi::op::sum<int>());
      log.Log0Verbose1() << program_timer.GetTimeString()
                         << " Done setting block id from Python function. "
                         << "Number of cells modified = " << global_num_cells_modified << ".";
    },
    R"(
    Set block ID from a function.

    Parameters
    ----------
    func: Callable[[pyopensn.math.Vector3, int], int]
        Function/lambda computing new block ID from cell centroid and old block ID.

    Examples
    --------
    >>> # Change block ID from 0 to 1 for cells inside the unit sphere.
    >>> def block_id_setter(cell_centroid, old_id):
    ...     if (old_id == 0) and (cell_centroid.Norm() < 1.0):
    ...         return 1
    ...     return old_id
    >>> mesh.SetBlockIDFromFunction(block_id_setter)
    )",
    py::arg("func")
  );
#ifdef __comment__  // to be removed once its uselessness is confirmed
  mesh_continuum.def(
    "SetBoundaryIDFromFunction",
    [](MeshContinuum& self, const std::function<std::string(Vector3, Vector3, std::uint64_t)>& func)
    {
      if (mpi_comm.size() != 1)
      {
        throw std::logic_error("This function can only be used in serial mode.");
      }
      log.Log0Verbose1() << program_timer.GetTimeString()
                         << " Setting boundary id from Python function.";
      int local_num_faces_modified = 0;
      // get boundary ID map
      std::map<uint64_t, std::string>& grid_boundary_id_map = self.GetBoundaryIDMap();
      // change local cells
      for (Cell& cell : self.local_cells)
      {
        for (CellFace& face : cell.faces)
        {
          if (!face.has_neighbor)
          {
            std::string boundary_name = func(face.centroid, face.normal, face.neighbor_id);
            std::uint64_t boundary_id = self.MakeBoundaryID(boundary_name);
            if (face.neighbor_id != boundary_id)
            {
              face.neighbor_id = boundary_id;
              ++local_num_faces_modified;
              if (grid_boundary_id_map.count(boundary_id) == 0)
              {
                grid_boundary_id_map[boundary_id] = boundary_name;
              }
            }
          }
        }
      }
      // change ghost cells
      std::vector<std::uint64_t> ghost_ids = self.cells.GetGhostGlobalIDs();
      for (std::uint64_t ghost_id : ghost_ids)
      {
        Cell& cell = self.cells[ghost_id];
        for (CellFace& face : cell.faces)
        {
          if (!face.has_neighbor)
          {
            std::string boundary_name = func(face.centroid, face.normal, face.neighbor_id);
            std::uint64_t boundary_id = self.MakeBoundaryID(boundary_name);
            if (face.neighbor_id != boundary_id)
            {
              face.neighbor_id = boundary_id;
              ++local_num_faces_modified;
              if (grid_boundary_id_map.count(boundary_id) == 0)
              {
                grid_boundary_id_map[boundary_id] = boundary_name;
              }
            }
          }
        }
      }
      // print number of modified cells
      int global_num_faces_modified;
      mpi_comm.all_reduce(local_num_faces_modified, global_num_faces_modified, mpi::op::sum<int>());
      log.Log0Verbose1() << program_timer.GetTimeString()
                         << " Done setting boundary id from lua function. "
                         << "Number of cells modified = " << global_num_faces_modified << ".";
    },
    R"(
    Set boundary ID from a function.

    Parameters
    ----------
    func: Callable[[pyopensn.math.Vector3, pyopensn.math.Vector3, int], str]
        Function/lambda computing new boundary name from face centroid, normal vector and current
        neighbor ID.
        (IDK how on earth users would know what boundary name corresponding to old neighbor ID)???
    )",
    py::arg("func")
  );
#endif  // __comment__

  // surface mesh
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
        Surface mesh filename.
    as_poly: bool, default=False
        Indicate if the surface mesh is allowed to contain polygonal facets (as opposed to only
        triangular faces).
    translation: pyopensn.math.Vector3, default=(0.0, 0.0, 0.0)
        Translation to perform on the mesh.
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
        Surface mesh filename.
    as_poly: bool
        Indicate if the surface mesh is allowed to contain polygonal facets (as opposed to only
        triangular faces).
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
        Surface mesh filename.
    as_poly: bool
        Indicate if the surface mesh is allowed to contain polygonal facets (as opposed to only
        triangular faces).
    )",
    py::arg("file_name"),
    py::arg("as_poly")
  );
  // clang-format on
}

// Wrap mesh generator
void
WrapMeshGenerator(py::module& mesh)
{
  // clang-format off
  // base mesh generator
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

  // extruded mesh generator
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
    R"(
    Construct an extruded mesh generator.

    Extrude 2D geometry using extrusion layers. Each layer is specified by either:

     - ``n`` and ``z`` to compute the z-levels automatically.
     - ``n`` and ``h`` to compute the h-levels automatically.

    The list of layers can be specified with a mixture of both ways.

    Parameters
    ----------
    scale: float, default=1.0
        Uniform scale to apply to the mesh after reading.
    inputs: List[pyopensn.mesh.MeshGenerator], default=[]
        A list of MeshGenerator objects.
    partitioner: List[pyopensn.mesh.GraphPartitioner], default=[]
        Handle to a GraphPartitioner object to use for parallel partitioning. This will default to
        PETScGraphPartitioner with a "parmetis" setting.
    replicated_mesh: bool, default=False
        Flag, when set, makes the mesh appear in full fidelity on each process.
    layers: List[Dict]
        List of layers. Parameters of each layers are represented as Python dictionary.
    top_boundary_name: str, default='ZMAX'
        The name to associate with the top boundary.
    bottom_boundary_name: str, default='ZMIN'
        The name to associate with the bottom boundary.

    Examples
    --------
    >>> emg = ExtruderMeshGenerator(
    ...     layers=[{"n": 1, "z": 2}, {"n": 2, "h": 0.5}]
    ... )
    )"
  );

  // orthogonal mesh generator
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
    R"(
    Construct an orthogonal mesh generator.

    Parameters
    ----------
    scale: float, default=1.0
        Uniform scale to apply to the mesh after reading.
    inputs: List[pyopensn.mesh.MeshGenerator], default=[]
        A list of MeshGenerator objects.
    partitioner: List[pyopensn.mesh.GraphPartitioner], default=[]
        Handle to a GraphPartitioner object to use for parallel partitioning. This will default to
        PETScGraphPartitioner with a "parmetis" setting.
    replicated_mesh: bool, default=False
        Flag, when set, makes the mesh appear in full fidelity on each process.
    node_sets: List[List[float]]
        Sets of nodes per dimension. Node values must be monotonically increasing.
    )"
  );

  // from file mesh generator
  auto from_file_mesh_generator = py::class_<FromFileMeshGenerator,
                                             std::shared_ptr<FromFileMeshGenerator>, MeshGenerator>(
    mesh,
    "FromFileMeshGenerator",
    R"(
    From file mesh generator.

    Wrapper of :cpp:class:`opensn::FromFileMeshGenerator`.
    )"
  );
  from_file_mesh_generator.def(
    py::init(
      [](py::kwargs & params)
      {
        return FromFileMeshGenerator::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a from-file mesh generator.

    Parameters
    ----------
    scale: float, default=1.0
        Uniform scale to apply to the mesh after reading.
    inputs: List[pyopensn.mesh.MeshGenerator], default=[]
        A list of MeshGenerator objects.
    partitioner: List[pyopensn.mesh.GraphPartitioner], default=[]
        Handle to a GraphPartitioner object to use for parallel partitioning. This will default to
        PETScGraphPartitioner with a "parmetis" setting.
    replicated_mesh: bool, default=False
        Flag, when set, makes the mesh appear in full fidelity on each process.
    filename: str
        Path to the file.
    block_id_fieldname: str, default='BlockID'
        The name of the field storing cell block/block ids. Only really used for .vtu, .pvtu and
        .e files.
    boundary_id_fieldname: str, default=''
        The name of the field storing boundary-ids.
    )"
  );

  // split file mesh generator
  auto split_file_mesh_generator = py::class_<SplitFileMeshGenerator,
                                              std::shared_ptr<SplitFileMeshGenerator>,
                                              MeshGenerator>(
    mesh,
    "SplitFileMeshGenerator",
    R"(
    Split file mesh generator.

    Wrapper of :cpp:class:`opensn::SplitFileMeshGenerator`.

    Generates the mesh only on location 0, thereafter partitions the mesh but instead of
    broadcasting the mesh to other locations it creates binary mesh files for each location.
    )"
  );
  split_file_mesh_generator.def(
    py::init(
      [](py::kwargs & params)
      {
        return SplitFileMeshGenerator::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a split-file mesh generator.

    Parameters
    ----------
    scale: float, default=1.0
        Uniform scale to apply to the mesh after reading.
    inputs: List[pyopensn.mesh.MeshGenerator], default=[]
        A list of MeshGenerator objects.
    partitioner: List[pyopensn.mesh.GraphPartitioner], default=[]
        Handle to a GraphPartitioner object to use for parallel partitioning. This will default to
        PETScGraphPartitioner with a "parmetis" setting.
    replicated_mesh: bool, default=False
        Flag, when set, makes the mesh appear in full fidelity on each process.
    num_partitions: int, default=0
        The number of partitions to generate. If zero, it will default to the number of MPI
        processes. Automatically ignored if the number of MPI processes greater 1.
    split_mesh_dir_path: str, default='split_mesh'
        Path of the directory to be created for containing the split meshes.
    file_prefix: str, default=''
        Prefix to use for all split mesh files. If not provided, it default to the input path's
        folder.
    read_only: bool, default=False
        Controls whether the split mesh is recreated or just read.
    verbosity_level: int, default=1
        Verbosity level. 1 will report each 10% complete. 2 will print each part and the number of
        local cells it wrote.
    )"
  );

  // distributed mesh generator
  auto distributed_mesh_generator = py::class_<DistributedMeshGenerator,
                                               std::shared_ptr<DistributedMeshGenerator>,
                                               MeshGenerator>(
    mesh,
    "DistributedMeshGenerator",
    R"(
    Distributed mesh generator.

    This class is responsible for generating a mesh, partitioning it, and distributing the
    individual partitions to different MPI locations. The mesh is generated on location 0,
    partitioned into multiple parts, serialized, and distributed to all other MPI ranks.

    Wrapper of :cpp:class:`opensn::DistributedMeshGenerator`.
    )"
  );
  distributed_mesh_generator.def(
    py::init(
      [](py::kwargs & params)
      {
        return DistributedMeshGenerator::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a distributed mesh generator.

    Parameters
    ----------
    scale: float, default=1.0
        Uniform scale to apply to the mesh after reading.
    inputs: List[pyopensn.mesh.MeshGenerator], default=[]
        A list of MeshGenerator objects.
    partitioner: List[pyopensn.mesh.GraphPartitioner], default=[]
        Handle to a GraphPartitioner object to use for parallel partitioning. This will default to
        PETScGraphPartitioner with a "parmetis" setting.
    replicated_mesh: bool, default=False
        Flag, when set, makes the mesh appear in full fidelity on each process.
    )"
  );
  // clang-format on
}

// Wrap graph partitioner
void
WrapGraphPartitioner(py::module& mesh)
{
  // clang-format off
  // base graph partitioner
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

  // linear graph partitioner
  auto linear_graph_partitioner = py::class_<LinearGraphPartitioner,
                                             std::shared_ptr<LinearGraphPartitioner>,
                                             GraphPartitioner>(
    mesh,
    "LinearGraphPartitioner",
    R"(
    Linear graph partitioner.

    Wrapper of :cpp:class:`opensn::LinearGraphPartitioner`.
    )"
  );
  linear_graph_partitioner.def(
    py::init(
      [](py::kwargs & params)
      {
        return LinearGraphPartitioner::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a linear graph partitioner.
    )"
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
  // clang-format on
}

void
py_mesh(py::module& pyopensn)
{
  py::module mesh = pyopensn.def_submodule("mesh", "Mesh generation module.");
  WrapMesh(mesh);
  WrapMeshGenerator(mesh);
  WrapGraphPartitioner(mesh);
}

} // namespace opensn
