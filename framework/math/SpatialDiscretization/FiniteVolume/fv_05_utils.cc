#include "FiniteVolume.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/UnknownManager/unknown_manager.h"

#include "chi_mpi.h"

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Get the number of ghost degrees-of-freedom.*/
size_t
chi_math::spatial_discretization::FiniteVolume::GetNumGhostDOFs(
  const chi_math::UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return ref_grid_.cells.GetNumGhosts() * N;
}

// ###################################################################
/**Returns the ghost DOF indices.*/
std::vector<int64_t>
chi_math::spatial_discretization::FiniteVolume::GetGhostDOFIndices(
  const chi_math::UnknownManager& unknown_manager) const
{
  std::vector<int64_t> dof_ids;
  dof_ids.reserve(GetNumGhostDOFs(unknown_manager));

  std::vector<uint64_t> ghost_cell_ids = ref_grid_.cells.GetGhostGlobalIDs();

  const size_t num_uks = unknown_manager.unknowns_.size();

  for (const auto cell_id : ghost_cell_ids)
  {
    const auto& cell = ref_grid_.cells[cell_id];
    for (size_t u = 0; u < num_uks; ++u)
    {
      const auto& unkn = unknown_manager.unknowns_[u];
      const size_t num_comps = unkn.num_components_;
      for (size_t c = 0; c < num_comps; ++c)
      {
        const int64_t dofmap = MapDOF(cell, 0, unknown_manager, u, c);
        dof_ids.push_back(dofmap);
      } // for c
    }   // for u
  }

  return dof_ids;
}

} // namespace chi_math::spatial_discretization
