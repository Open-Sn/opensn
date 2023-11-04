#include "MeshGenerator.h"

#include "mesh/UnpartitionedMesh/unpartitioned_mesh.h"
#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/VolumeMesher/chi_volumemesher.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "graphs/GraphPartitioner.h"
#include "graphs/PETScGraphPartitioner.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_mesh
{

RegisterChiObject(chi_mesh, MeshGenerator);

chi::InputParameters
MeshGenerator::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription("The base class for all mesh generators");
  params.SetDocGroup("doc_MeshGenerators");

  params.AddOptionalParameter("scale", 1.0, "Uniform scale to apply to the mesh after reading.");

  params.AddOptionalParameterArray(
    "inputs", std::vector<size_t>{}, "A list of handles to MeshGenerator objects");

  params.AddOptionalParameter(
    "partitioner",
    0,
    "Handle to a GraphPartitioner object to use for parallel partitioning."
    "This will default to PETScGraphPartitioner with a \"parmetis\" setting");

  params.AddOptionalParameter(
    "replicated_mesh",
    false,
    "Flag, when set, makes the mesh appear in full fidelity on each process");

  return params;
}

MeshGenerator::MeshGenerator(const chi::InputParameters& params)
  : ChiObject(params),
    scale_(params.GetParamValue<double>("scale")),
    replicated_(params.GetParamValue<bool>("replicated_mesh"))
{
  //============================================= Convert input handles
  auto input_handles = params.GetParamVectorValue<size_t>("inputs");

  for (const size_t input_handle : input_handles)
  {
    auto& mesh_generator =
      Chi::GetStackItem<MeshGenerator>(Chi::object_stack, input_handle, __FUNCTION__);
    inputs_.push_back(&mesh_generator);
  }

  //============================================= Set partitioner
  size_t partitioner_handle;
  if (params.ParametersAtAssignment().Has("partitioner"))
    partitioner_handle = params.GetParamValue<size_t>("partitioner");
  else
  {
    auto& factory = ChiObjectFactory::GetInstance();
    auto valid_params = chi::PETScGraphPartitioner::GetInputParameters();
    partitioner_handle =
      factory.MakeRegisteredObjectOfType("chi::PETScGraphPartitioner", chi::ParameterBlock());
  }
  partitioner_ =
    &Chi::GetStackItem<chi::GraphPartitioner>(Chi::object_stack, partitioner_handle, __FUNCTION__);
}

// ##################################################################
/**Default behavior here is to return the input umesh unaltered.*/
std::unique_ptr<UnpartitionedMesh>
MeshGenerator::GenerateUnpartitionedMesh(std::unique_ptr<UnpartitionedMesh> input_umesh)
{
  return input_umesh;
}

/**Final execution step. */
void
MeshGenerator::Execute()
{
  //======================================== Execute all input generators
  // Note these could be empty
  std::unique_ptr<UnpartitionedMesh> current_umesh = nullptr;
  for (auto mesh_generator_ptr : inputs_)
  {
    auto new_umesh = mesh_generator_ptr->GenerateUnpartitionedMesh(std::move(current_umesh));
    current_umesh = std::move(new_umesh);
  }

  //======================================== Generate final umesh and convert it
  current_umesh = GenerateUnpartitionedMesh(std::move(current_umesh));

  std::vector<int64_t> cell_pids;
  if (Chi::mpi.location_id == 0) cell_pids = PartitionMesh(*current_umesh, Chi::mpi.process_count);

  BroadcastPIDs(cell_pids, 0, Chi::mpi.comm);

  auto grid_ptr = SetupMesh(std::move(current_umesh), cell_pids);

  //======================================== Assign the mesh to a VolumeMesher
  auto new_mesher = std::make_shared<chi_mesh::VolumeMesher>(VolumeMesherType::UNPARTITIONED);
  new_mesher->SetContinuum(grid_ptr);

  if (Chi::current_mesh_handler < 0) chi_mesh::PushNewHandlerAndGetIndex();

  auto& cur_hndlr = chi_mesh::GetCurrentHandler();
  cur_hndlr.SetVolumeMesher(new_mesher);

  Chi::mpi.Barrier();
}

void
MeshGenerator::SetGridAttributes(chi_mesh::MeshContinuum& grid,
                                 MeshAttributes new_attribs,
                                 std::array<size_t, 3> ortho_cells_per_dimension)
{
  grid.SetAttributes(new_attribs, ortho_cells_per_dimension);
}

void
MeshGenerator::ComputeAndPrintStats(const chi_mesh::MeshContinuum& grid)
{
  const size_t num_local_cells = grid.local_cells.size();
  size_t num_global_cells = 0;

  MPI_Allreduce(&num_local_cells,       // sendbuf
                &num_global_cells,      // recvbuf
                1,                      // count
                MPI_UNSIGNED_LONG_LONG, // datatype
                MPI_SUM,                // operation
                Chi::mpi.comm);         // communicator

  size_t max_num_local_cells;
  MPI_Allreduce(&num_local_cells,       // sendbuf
                &max_num_local_cells,   // recvbuf
                1,                      // count
                MPI_UNSIGNED_LONG_LONG, // datatype
                MPI_MAX,                // operation
                Chi::mpi.comm);         // communicator

  size_t min_num_local_cells;
  MPI_Allreduce(&num_local_cells,       // sendbuf
                &min_num_local_cells,   // recvbuf
                1,                      // count
                MPI_UNSIGNED_LONG_LONG, // datatype
                MPI_MIN,                // operation
                Chi::mpi.comm);         // communicator

  const size_t avg_num_local_cells = num_global_cells / Chi::mpi.process_count;
  const size_t num_local_ghosts = grid.cells.GetNumGhosts();
  const double local_ghost_to_local_cell_ratio = double(num_local_ghosts) / double(num_local_cells);

  double average_ghost_ratio;
  MPI_Allreduce(&local_ghost_to_local_cell_ratio, // sendbuf
                &average_ghost_ratio,             // recvbuf
                1,                                // count
                MPI_DOUBLE,                       // datatype
                MPI_SUM,                          // operation
                Chi::mpi.comm);                   // communicator

  average_ghost_ratio /= Chi::mpi.process_count;

  std::stringstream outstr;
  outstr << "Mesh statistics:\n";
  outstr << "  Global cell count             : " << num_global_cells << "\n";
  outstr << "  Local cell count (avg,max,min): ";
  outstr << avg_num_local_cells << ",";
  outstr << max_num_local_cells << ",";
  outstr << min_num_local_cells << "\n";
  outstr << "  Ghost-to-local ratio (avg)    : " << average_ghost_ratio;

  Chi::log.Log() << "\n" << outstr.str() << "\n\n";

  Chi::log.LogAllVerbose2() << Chi::mpi.location_id << "Local cells=" << num_local_cells;
}

} // namespace chi_mesh
