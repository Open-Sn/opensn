#include "PETScGraphPartitioner.h"

#include "ChiObjectFactory.h"

#include "petsc.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi
{

RegisterChiObject(chi, PETScGraphPartitioner);

InputParameters PETScGraphPartitioner::GetInputParameters()
{
  InputParameters params = GraphPartitioner::GetInputParameters();

  params.SetGeneralDescription("PETSc based partitioning");
  params.SetDocGroup("Graphs");

  params.AddOptionalParameter(
    "type", "parmetis", "The type of PETSc partitioner");

  return params;
}

PETScGraphPartitioner::PETScGraphPartitioner(const InputParameters& params)
  : GraphPartitioner(params), type_(params.GetParamValue<std::string>("type"))
{
}

std::vector<int64_t> PETScGraphPartitioner::Partition(
  const std::vector<std::vector<uint64_t>>& graph,
  const std::vector<chi_mesh::Vector3>&,
  int number_of_parts)
{
  Chi::log.Log0Verbose1() << "Partitioning with PETScGraphPartitioner";
  //================================================== Determine avg num faces
  //                                                   per cell
  // This is done so we can reserve size better
  const size_t num_raw_cells = graph.size();
  size_t num_raw_faces = 0;
  for (auto& cell_row : graph)
    num_raw_faces += cell_row.size();
  size_t avg_num_face_per_cell = std::ceil(static_cast<double>(num_raw_faces) /
                                           static_cast<double>(num_raw_cells));

  //================================================== Start building indices
  std::vector<int64_t> cell_pids(num_raw_cells, 0);
  if (num_raw_cells > 1)
  {
    //======================================== Build indices
    std::vector<int64_t> i_indices(num_raw_cells + 1, 0);
    std::vector<int64_t> j_indices;
    j_indices.reserve(num_raw_cells * avg_num_face_per_cell);
    {
      int64_t i = 0;
      int64_t icount = 0;
      for (const auto& cell : graph)
      {
        i_indices[i] = icount;

        for (const uint64_t neighbor_id : cell)
        {
          j_indices.push_back(static_cast<int64_t>(neighbor_id));
          ++icount;
        }
        ++i;
      }
      i_indices[i] = icount;
    }

    Chi::log.Log0Verbose1() << "Done building indices.";

    //======================================== Copy to raw arrays
    int64_t* i_indices_raw;
    int64_t* j_indices_raw;
    PetscMalloc(i_indices.size() * sizeof(int64_t), &i_indices_raw);
    PetscMalloc(j_indices.size() * sizeof(int64_t), &j_indices_raw);

    for (int64_t j = 0; j < static_cast<int64_t>(i_indices.size()); ++j)
      i_indices_raw[j] = i_indices[j];

    for (int64_t j = 0; j < static_cast<int64_t>(j_indices.size()); ++j)
      j_indices_raw[j] = j_indices[j];

    Chi::log.Log0Verbose1() << "Done copying to raw indices.";

    //========================================= Create adjacency matrix
    Mat Adj; // Adjacency matrix
    MatCreateMPIAdj(PETSC_COMM_SELF,
                    (int64_t)num_raw_cells,
                    (int64_t)num_raw_cells,
                    i_indices_raw,
                    j_indices_raw,
                    nullptr,
                    &Adj);

    Chi::log.Log0Verbose1() << "Done creating adjacency matrix.";

    //========================================= Create partitioning
    MatPartitioning part;
    IS is, isg;
    MatPartitioningCreate(MPI_COMM_SELF, &part);
    MatPartitioningSetAdjacency(part, Adj);
    MatPartitioningSetType(part, type_.c_str());
    MatPartitioningSetNParts(part, number_of_parts);
    MatPartitioningApply(part, &is);
    MatPartitioningDestroy(&part);
    MatDestroy(&Adj);
    ISPartitioningToNumbering(is, &isg);
    Chi::log.Log0Verbose1() << "Done building paritioned index set.";

    //========================================= Get cell global indices
    const int64_t* cell_pids_raw;
    ISGetIndices(is, &cell_pids_raw);
    for (size_t i = 0; i < num_raw_cells; ++i)
      cell_pids[i] = cell_pids_raw[i];
    ISRestoreIndices(is, &cell_pids_raw);

    Chi::log.Log0Verbose1() << "Done retrieving cell global indices.";
  } // if more than 1 cell

  Chi::log.Log0Verbose1() << "Done partitioning with PETScGraphPartitioner";
  return cell_pids;
}

} // namespace chi