#include "FiniteVolume.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/UnknownManager/unknown_manager.h"

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Builds finite volume based sparsity pattern.*/
void
FiniteVolume::BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                                   std::vector<int64_t>& nodal_nnz_off_diag,
                                   const chi_math::UnknownManager& unknown_manager) const
{
  unsigned int num_uk = unknown_manager.unknowns_.size();          // Number of unknowns
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize(); // Total number of unknowns

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  const size_t num_local_cells = ref_grid_.local_cells.size();

  nodal_nnz_in_diag.resize(num_local_cells * N, 0.0);
  nodal_nnz_off_diag.resize(num_local_cells * N, 0.0);

  for (int uk = 0; uk < num_uk; ++uk)
  {
    const unsigned int num_comps = unknown_manager.unknowns_[uk].num_components_;
    for (int comp = 0; comp < num_comps; ++comp)
    {
      for (auto& cell : ref_grid_.local_cells)
      {
        const int64_t i = MapDOFLocal(cell, 0, unknown_manager, uk, comp);

        nodal_nnz_in_diag[i] += 1;

        for (auto& face : cell.faces_)
        {
          if (not face.has_neighbor_) continue;

          if (face.IsNeighborLocal(ref_grid_)) nodal_nnz_in_diag[i] += 1;
          else
            nodal_nnz_off_diag[i] += 1;
        }
      } // for cell
    }   // for components
  }     // for unknown
}

} // namespace chi_math::spatial_discretization
