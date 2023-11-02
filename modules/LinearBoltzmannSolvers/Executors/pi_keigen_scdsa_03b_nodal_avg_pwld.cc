#include "pi_keigen_scdsa.h"

#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

namespace lbs
{

/**This method takes an input vector that is the local version of
 * a PWLD discrete space and then makes it continuous by applying nodal
 * averages.*/
std::vector<double>
XXPowerIterationKEigenSCDSA::NodallyAveragedPWLDVector(
  const std::vector<double>& input,
  const chi_math::SpatialDiscretization& pwld_sdm,
  const chi_math::SpatialDiscretization& pwlc_sdm,
  const chi_math::UnknownManager& uk_man,
  const XXPowerIterationKEigenSCDSA::GhostInfo& ghost_info)
{
  const auto& vgc = ghost_info.vector_ghost_communicator;
  const auto& dfem_dof_global2local_map = ghost_info.ghost_global_id_2_local_map;

  auto input_with_ghosts = vgc->MakeGhostedVector(input);
  vgc->CommunicateGhostEntries(input_with_ghosts);

  typedef const int64_t cint64_t;

  const auto& grid = pwld_sdm.Grid();

  const size_t num_unknowns = uk_man.unknowns_.size();

  const size_t num_cfem_local_dofs = pwlc_sdm.GetNumLocalAndGhostDOFs(uk_man);

  std::vector<double> cont_input(num_cfem_local_dofs, 0.0);
  std::vector<double> cont_input_ctr(num_cfem_local_dofs, 0.0);

  std::map<int64_t, int64_t> cfem_dof_global2local_map;

  //================================================== Local cells first
  std::set<uint64_t> partition_bndry_vertex_id_set;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns_[u].num_components_;
        for (size_t c = 0; c < num_components; ++c)
        {
          cint64_t dof_dfem_map = pwld_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          cint64_t dof_cfem_map = pwlc_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          cint64_t dof_cfem_map_globl = pwlc_sdm.MapDOF(cell, i, uk_man, u, c);

          cfem_dof_global2local_map[dof_cfem_map_globl] = dof_cfem_map;

          const double phi_value = input[dof_dfem_map];

          cont_input[dof_cfem_map] += phi_value;
          cont_input_ctr[dof_cfem_map] += 1.0;
        } // for component c
      }   // for unknown u
    }     // for node i

    for (const auto& face : cell.faces_)
      if (face.has_neighbor_)
        if (not grid.IsCellLocal(face.neighbor_id_))
          for (const uint64_t vid : face.vertex_ids_)
            partition_bndry_vertex_id_set.insert(vid);
  } // for local cell

  //================================================== Ghost cells
  const auto ghost_cell_ids = grid.cells.GetGhostGlobalIDs();
  const auto& vid_set = partition_bndry_vertex_id_set;
  for (const auto global_id : ghost_cell_ids)
  {
    const auto& cell = grid.cells[global_id];
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      if (vid_set.find(cell.vertex_ids_[i]) == vid_set.end()) continue;

      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns_[u].num_components_;
        for (size_t c = 0; c < num_components; ++c)
        {
          cint64_t dof_dfem_map_globl = pwld_sdm.MapDOF(cell, i, uk_man, u, c);
          cint64_t dof_cfem_map_globl = pwlc_sdm.MapDOF(cell, i, uk_man, u, c);
          if (cfem_dof_global2local_map.count(dof_cfem_map_globl) > 0)
          {
            cint64_t dof_dfem_map = dfem_dof_global2local_map.at(dof_dfem_map_globl);
            cint64_t dof_cfem_map = cfem_dof_global2local_map[dof_cfem_map_globl];

            const double phi_value = input_with_ghosts[dof_dfem_map];

            cont_input[dof_cfem_map] += phi_value;
            cont_input_ctr[dof_cfem_map] += 1.0;
          }
        } // for component
      }   // for unknown
    }     // for node i
  }       // for ghost cell

  //================================================== Compute nodal averages
  {
    const size_t num_vals = cont_input.size();
    for (size_t k = 0; k < num_vals; ++k)
      cont_input[k] /= cont_input_ctr[k];
  }

  //================================================== Project back to dfem
  auto output = input;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns_[u].num_components_;
        for (size_t c = 0; c < num_components; ++c)
        {
          cint64_t dof_dfem_map = pwld_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          cint64_t dof_cfem_map = pwlc_sdm.MapDOFLocal(cell, i, uk_man, u, c);

          const double phi_value = cont_input[dof_cfem_map];

          output[dof_dfem_map] = phi_value;
        } // for component c
      }   // for unknown u
    }     // for node i
  }       // for local cell

  return output;
}

} // namespace lbs