#include "volmesher_predefunpart.h"

#include "chi_runtime.h"
#include "chi_mpi.h"


//###################################################################
/**Determines if a chi_mesh::UnpartitionedMesh::LightWeightCell is a
 * neighbor to the current partition for ParMETIS-style partitioning.
 * This method loops over the faces of the lightweight cell and
 * determines the partition-id of each the neighbors. If the neighbor
 * has a partition id equal to that of the current process then
 * it means this reference cell is a neighbor.*/
bool chi_mesh::VolumeMesherPredefinedUnpartitioned::
  CellHasLocalScope(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
    uint64_t cell_global_id,
    const std::vector<std::set<uint64_t>>& vertex_subscriptions,
    const std::vector<int64_t>& cell_partition_ids)
{
  //First determine if the cell is a local cell
  int cell_pid = static_cast<int>(cell_partition_ids[cell_global_id]);
  if (cell_pid == Chi::mpi.location_id)
    return true;

  //Now determine if the cell is a ghost cell
  for (uint64_t vid : lwcell.vertex_ids)
    for (uint64_t cid : vertex_subscriptions[vid])
    {
      if (cid == cell_global_id) continue;
      int adj_pid = static_cast<int>(cell_partition_ids[cid]);
      if (adj_pid == Chi::mpi.location_id)
        return true;
    }

  return false;
}