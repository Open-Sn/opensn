#include "volmesher_predefunpart.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"


#include "petsc.h"

//###################################################################
/** Applies KBA-style partitioning to the mesh.*/
std::vector<int64_t> chi_mesh::VolumeMesherPredefinedUnpartitioned::
  PARMETIS(const UnpartitionedMesh &umesh)
{
  Chi::log.Log() << "Partitioning mesh with ParMETIS.";

  //================================================== Determine avg num faces
  //                                                   per cell
  const size_t num_raw_cells = umesh.GetNumberOfCells();
  size_t num_raw_faces = 0;
  for (auto& cell : umesh.GetRawCells())
    num_raw_faces += cell->faces.size();
  size_t avg_num_face_per_cell =
    std::ceil(static_cast<double>(num_raw_faces)/
              static_cast<double>(num_raw_cells));

  //================================================== Start building indices
  std::vector<int64_t> cell_pids(num_raw_cells, 0);
  if (Chi::mpi.location_id == 0)
  {
    if (num_raw_cells > 1)
    {
      //======================================== Build indices
      std::vector<int64_t> i_indices(num_raw_cells+1,0);
      std::vector<int64_t> j_indices;
      j_indices.reserve(num_raw_cells * avg_num_face_per_cell);
      {
        int64_t i=0;
        int64_t icount = 0;
        for (auto cell : umesh.GetRawCells())
        {
          i_indices[i] = icount;

          for (auto& face : cell->faces)
            if (face.has_neighbor)
            {
              j_indices.push_back(static_cast<int64_t>(face.neighbor));
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
      PetscMalloc(i_indices.size()*sizeof(int64_t),&i_indices_raw);
      PetscMalloc(j_indices.size()*sizeof(int64_t),&j_indices_raw);

      for (int64_t j=0; j<static_cast<int64_t>(i_indices.size()); ++j)
        i_indices_raw[j] = i_indices[j];

      for (int64_t j=0; j<static_cast<int64_t>(j_indices.size()); ++j)
        j_indices_raw[j] = j_indices[j];

      Chi::log.Log0Verbose1() << "Done copying to raw indices.";

      //========================================= Create adjacency matrix
      Mat Adj; //Adjacency matrix
      MatCreateMPIAdj(PETSC_COMM_SELF,
                      (int64_t)num_raw_cells,
                      (int64_t)num_raw_cells,
                      i_indices_raw, j_indices_raw, nullptr, &Adj);

      Chi::log.Log0Verbose1() << "Done creating adjacency matrix.";

      //========================================= Create partitioning
      MatPartitioning part;
      IS is,isg;
      MatPartitioningCreate(MPI_COMM_SELF,&part);
      MatPartitioningSetAdjacency(part,Adj);
      MatPartitioningSetType(part,"parmetis");
      MatPartitioningSetNParts(part, Chi::mpi.process_count);
      MatPartitioningApply(part,&is);
      MatPartitioningDestroy(&part);
      MatDestroy(&Adj);
      ISPartitioningToNumbering(is,&isg);
      Chi::log.Log0Verbose1() << "Done building paritioned index set.";

      //========================================= Get cell global indices
      const int64_t* cell_pids_raw;
      ISGetIndices(is,&cell_pids_raw);
      for (size_t i=0; i<num_raw_cells; ++i)
        cell_pids[i] = cell_pids_raw[i];
      ISRestoreIndices(is,&cell_pids_raw);

      Chi::log.Log0Verbose1() << "Done retrieving cell global indices.";
    }//if more than 1 cell
  }//if home location

  //======================================== Broadcast partitioning to all
  //                                         locations
  MPI_Bcast(cell_pids.data(),                 //buffer [IN/OUT]
            static_cast<int>(num_raw_cells),  //count
            MPI_LONG_LONG_INT,                //data type
            0,                                //root
            Chi::mpi.comm);                  //communicator
  Chi::log.Log() << "Done partitioning mesh.";

  return cell_pids;
}