#include "pi_keigen_scdsa.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

namespace lbs
{

// ##################################################################
/**Copies only the scalar moments from an lbs primary flux moments
 * vector.*/
std::vector<double>
XXPowerIterationKEigenSCDSA::CopyOnlyPhi0(const LBSGroupset& groupset,
                                          const std::vector<double>& phi_in)
{
  typedef const int64_t cint64;

  const auto& lbs_sdm = lbs_solver_.SpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();
  const auto& phi_uk_man = lbs_solver_.UnknownManager();

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  const size_t diff_num_local_dofs =
    requires_ghosts_ ? diff_sdm.GetNumLocalAndGhostDOFs(diff_uk_man)
                     : diff_sdm.GetNumLocalDOFs(diff_uk_man);

  std::vector<double> phi_data;
  if (continuous_sdm_ptr_)
    phi_data = NodallyAveragedPWLDVector(
      phi_in, lbs_sdm, diff_sdm, phi_uk_man, lbs_pwld_ghost_info_);
  else
    phi_data = phi_in;

  VecDbl output_phi_local(diff_num_local_dofs, 0.0);

  for (const auto& cell : lbs_solver_.Grid().local_cells)
  {
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      cint64 diff_phi_map = diff_sdm.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      cint64 lbs_phi_map = lbs_sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* output_mapped = &output_phi_local[diff_phi_map];
      const double* phi_in_mapped = &phi_data[lbs_phi_map];

      for (size_t g = 0; g < gss; g++)
      {
        output_mapped[g] = phi_in_mapped[g];
      } // for g
    }   // for node
  }     // for cell

  return output_phi_local;
}

// ##################################################################
/**Copies back only the scalar moments to a lbs primary flux vector.*/
void XXPowerIterationKEigenSCDSA::ProjectBackPhi0(
  const LBSGroupset& groupset,
  const std::vector<double>& input,
  std::vector<double>& output)
{
  typedef const int64_t cint64;

  const auto& lbs_sdm = lbs_solver_.SpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();
  const auto& phi_uk_man = lbs_solver_.UnknownManager();

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  const size_t diff_num_local_dofs =
    requires_ghosts_ ? diff_sdm.GetNumLocalAndGhostDOFs(diff_uk_man)
                     : diff_sdm.GetNumLocalDOFs(diff_uk_man);

  ChiLogicalErrorIf(input.size() != diff_num_local_dofs,
                    "Vector size mismatch");

  for (const auto& cell : lbs_solver_.Grid().local_cells)
  {
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      cint64 diff_phi_map = diff_sdm.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      cint64 lbs_phi_map = lbs_sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* input_mapped = &input[diff_phi_map];
      double* output_mapped = &output[lbs_phi_map];

      for (int g = 0; g < gss; g++)
        output_mapped[g] = input_mapped[g];
    } // for dof
  }   // for cell
}

// ##################################################################
/**Creates a ghost communicator and all associated information.*/
XXPowerIterationKEigenSCDSA::GhostInfo
XXPowerIterationKEigenSCDSA::MakePWLDVecGhostCommInfo(
  const chi_math::SpatialDiscretization& sdm,
  const chi_math::UnknownManager& uk_man)
{
  Chi::log.Log() << "Making PWLD ghost communicator";

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(uk_man);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(uk_man);

  Chi::log.Log() << "Number of global dofs" << num_globl_dofs;

  const size_t num_unknowns = uk_man.unknowns_.size();

  // Build a list of global ids
  std::set<int64_t> global_dof_ids_set;

  const auto& grid = lbs_solver_.Grid();
  const auto ghost_cell_ids = grid.cells.GetGhostGlobalIDs();
  for (const auto global_id : ghost_cell_ids)
  {
    const auto& cell = grid.cells[global_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_comps = uk_man.unknowns_[u].num_components_;
        for (size_t c = 0; c < num_comps; ++c)
        {
          const int64_t dof_map = sdm.MapDOF(cell, i, uk_man, u, c);
          global_dof_ids_set.insert(dof_map);
        } // for component
      }   // for unknown
    }     // for node i
  }       // for ghost cell

  // Convert the list to a vector
  std::vector<int64_t> global_indices(global_dof_ids_set.begin(),
                                      global_dof_ids_set.end());

  // Create the vector ghost communicator
  auto vgc = std::make_shared<chi_math::VectorGhostCommunicator>(
    num_local_dofs, num_globl_dofs, global_indices, Chi::mpi.comm);

  // Create the map
  std::map<int64_t, int64_t> ghost_global_id_2_local_map;
  {
    int64_t k = 0;
    for (const auto ghost_id : global_indices)
    {
      ghost_global_id_2_local_map[ghost_id] =
        static_cast<int64_t>(num_local_dofs + k++);
    }
  }

  Chi::log.Log() << "Done making PWLD ghost communicator";
  return {vgc, ghost_global_id_2_local_map};
}

} // namespace lbs