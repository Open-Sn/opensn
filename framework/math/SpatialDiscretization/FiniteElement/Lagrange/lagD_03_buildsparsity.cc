#include "LagrangeDiscontinuous.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Builds the sparsity pattern for a Discontinuous Finite Element Method.*/
void
LagrangeDiscontinuous::BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                                            std::vector<int64_t>& nodal_nnz_off_diag,
                                            const chi_math::UnknownManager& unknown_manager) const
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL CONNECTIVITY
  size_t local_dof_count = local_base_block_size_;

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag.resize(local_dof_count, 0);
  nodal_nnz_off_diag.resize(local_dof_count, 0);

  int lc = 0;
  for (const auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);
    size_t num_nodes = cell_mapping.NumNodes();

    //==================================== Self connection
    for (int i = 0; i < num_nodes; ++i)
    {
      int64_t ir = cell_local_block_address_[lc] + i;
      nodal_nnz_in_diag[ir] += static_cast<int64_t>(num_nodes);
    }

    //==================================== Local adjacent cell connections
    for (auto& face : cell.faces_)
    {
      if (face.has_neighbor_ and face.IsNeighborLocal(ref_grid_))
      {
        const auto& adj_cell = ref_grid_.cells[face.neighbor_id_];
        const auto& adj_cell_mapping = GetCellMapping(adj_cell);

        for (int i = 0; i < num_nodes; ++i)
        {
          int64_t ir = cell_local_block_address_[lc] + i;
          nodal_nnz_in_diag[ir] += static_cast<int64_t>(adj_cell_mapping.NumNodes());
        }
      }
    } // for face
    ++lc;
  } // for local cell

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORING CONNECTIVITY
  lc = 0;
  for (auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);

    //==================================== Local adjacent cell connections
    for (auto& face : cell.faces_)
    {
      if (face.has_neighbor_ and (not face.IsNeighborLocal(ref_grid_)))
      {
        const auto& adj_cell = ref_grid_.cells[face.neighbor_id_];
        const auto& adj_cell_mapping = GetCellMapping(adj_cell);

        for (int i = 0; i < cell_mapping.NumNodes(); ++i)
        {
          int64_t ir = cell_local_block_address_[lc] + i;
          nodal_nnz_off_diag[ir] += static_cast<int64_t>(adj_cell_mapping.NumNodes());
        }
      }
    }
    ++lc;
  } // for local cell

  //============================================= Spacing according to unknown
  //                                              manager
  auto backup_nnz_in_diag = nodal_nnz_in_diag;
  auto backup_nnz_off_diag = nodal_nnz_off_diag;

  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag.resize(local_base_block_size_ * N, 0);
  nodal_nnz_off_diag.resize(local_base_block_size_ * N, 0);

  if (unknown_manager.dof_storage_type_ == chi_math::UnknownStorageType::NODAL)
  {
    int ir = -1;
    for (int i = 0; i < local_base_block_size_; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        ++ir;
        nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
        nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
      } // for j
    }   // for i
  }
  else if (unknown_manager.dof_storage_type_ == chi_math::UnknownStorageType::BLOCK)
  {
    int ir = -1;
    for (int j = 0; j < N; ++j)
    {
      for (int i = 0; i < local_base_block_size_; ++i)
      {
        ++ir;
        nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
        nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
      } // for i
    }   // for j
  }

  Chi::mpi.Barrier();
}

} // namespace chi_math::spatial_discretization