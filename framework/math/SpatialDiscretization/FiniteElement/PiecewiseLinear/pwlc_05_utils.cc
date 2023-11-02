#include "PieceWiseLinearContinuous.h"

#include "chi_runtime.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#define sc_int64 static_cast<int64_t>

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Get the number of ghost degrees-of-freedom.*/
size_t PieceWiseLinearContinuous::GetNumGhostDOFs(
  const chi_math::UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return ghost_node_mapping_.size() * N;
}

// ###################################################################
/**Returns the ghost DOF indices.*/
std::vector<int64_t> PieceWiseLinearContinuous::GetGhostDOFIndices(
  const chi_math::UnknownManager& unknown_manager) const
{
  std::vector<int64_t> dof_ids;
  dof_ids.reserve(GetNumGhostDOFs(unknown_manager));

  const size_t num_unknown_comps =
    unknown_manager.GetTotalUnknownStructureSize();
  const auto storage = unknown_manager.dof_storage_type_;
  const size_t num_unknowns = unknown_manager.unknowns_.size();

  for (const auto& vid_gnid : ghost_node_mapping_)
  {
    const int64_t global_id = vid_gnid.second;

    for (size_t u = 0; u < num_unknowns; ++u)
    {
      const auto& unkn = unknown_manager.unknowns_[u];
      const size_t num_comps = unkn.num_components_;
      for (size_t c = 0; c < num_comps; ++c)
      {
        size_t block_id = unknown_manager.MapUnknown(u, c);
        int64_t address = -1;
        if (storage == chi_math::UnknownStorageType::BLOCK)
        {
          for (int locJ = 0; locJ < Chi::mpi.process_count; ++locJ)
          {
            const int64_t local_id =
              global_id - sc_int64(locJ_block_address_[locJ]);

            if (local_id < 0 or local_id >= locJ_block_size_[locJ]) continue;

            address = sc_int64(locJ_block_address_[locJ] * num_unknown_comps) +
                      sc_int64(locJ_block_size_[locJ] * block_id) + local_id;
            break;
          }
        }
        else if (storage == chi_math::UnknownStorageType::NODAL)
          address =
            global_id * sc_int64(num_unknown_comps) + sc_int64(block_id);

        dof_ids.push_back(address);
      } // for c
    }   // for u
  }

  return dof_ids;
}

} // namespace chi_math::spatial_discretization
