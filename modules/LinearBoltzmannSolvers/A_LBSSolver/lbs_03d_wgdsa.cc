#include "A_LBSSolver/lbs_solver.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "A_LBSSolver/Acceleration/diffusion_mip_solver.h"

// ###################################################################
/**Initializes the Within-Group DSA solver. */
void
lbs::LBSSolver::InitWGDSA(LBSGroupset& groupset, bool vaccum_bcs_are_dirichlet /*=true*/)
{
  if (groupset.apply_wgdsa_)
  {
    //=========================================== Make UnknownManager
    const size_t num_gs_groups = groupset.groups_.size();
    chi_math::UnknownManager uk_man;
    uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_gs_groups);

    //=========================================== Make boundary conditions
    auto bcs = acceleration::TranslateBCs(sweep_boundaries_, vaccum_bcs_are_dirichlet);

    //=========================================== Make xs map
    auto matid_2_mgxs_map = acceleration::PackGroupsetXS(
      matid_to_xs_map_, groupset.groups_.front().id_, groupset.groups_.back().id_);

    //=========================================== Create solver
    const auto& sdm = *discretization_;

    auto solver =
      std::make_shared<acceleration::DiffusionMIPSolver>(std::string(TextName() + "_WGDSA"),
                                                         sdm,
                                                         uk_man,
                                                         bcs,
                                                         matid_2_mgxs_map,
                                                         unit_cell_matrices_,
                                                         true); // verbosity
    chi::ParameterBlock block;

    solver->options.residual_tolerance = groupset.wgdsa_tol_;
    solver->options.max_iters = groupset.wgdsa_max_iters_;
    solver->options.verbose = groupset.wgdsa_verbose_;
    solver->options.additional_options_string = groupset.wgdsa_string_;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);

    solver->AssembleAand_b(dummy_rhs);

    groupset.wgdsa_solver_ = solver;
  }
}

// ###################################################################
/**Cleans up memory consuming items. */
void
lbs::LBSSolver::CleanUpWGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_wgdsa_) groupset.wgdsa_solver_ = nullptr;
}

// ###################################################################
/**Creates a vector from a lbs primary stl vector where only the
 * scalar moments are mapped to the DOFs needed by WGDSA.*/
std::vector<double>
lbs::LBSSolver::WGSCopyOnlyPhi0(const LBSGroupset& groupset, const std::vector<double>& phi_in)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  std::vector<double> output_phi_local(sdm.GetNumLocalDOFs(dphi_uk_man), 0.0);

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* output_mapped = &output_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g = 0; g < gss; g++)
      {
        output_mapped[g] = phi_in_mapped[g];
      } // for g
    }   // for node
  }     // for cell

  return output_phi_local;
}

// ###################################################################
/**From the WGDSA DOFs, projects the scalar moments back into a
 * primary STL vector.*/
void
lbs::LBSSolver::GSProjectBackPhi0(const LBSGroupset& groupset,
                                  const std::vector<double>& input,
                                  std::vector<double>& output)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* input_mapped = &input[dphi_map];
      double* output_mapped = &output[phi_map];

      for (int g = 0; g < gss; g++)
        output_mapped[g] = input_mapped[g];
    } // for dof
  }   // for cell
}

// ###################################################################
/**Assembles a delta-phi vector on the first moment.*/
void
lbs::LBSSolver::AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                            const std::vector<double>& phi_in,
                                            std::vector<double>& delta_phi_local)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  delta_phi_local.clear();
  delta_phi_local.assign(sdm.GetNumLocalDOFs(dphi_uk_man), 0.0);

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& sigma_s = matid_to_xs_map_[cell.material_id_]->SigmaSGtoG();

    for (size_t i = 0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* delta_phi_mapped = &delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g = 0; g < gss; g++)
      {
        delta_phi_mapped[g] = sigma_s[gsi + g] * phi_in_mapped[g];
      } // for g
    }   // for node
  }     // for cell
}

// ###################################################################
/**DAssembles a delta-phi vector on the first moment.*/
void
lbs::LBSSolver::DisAssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                               const std::vector<double>& delta_phi_local,
                                               std::vector<double>& ref_phi_new)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* delta_phi_mapped = &delta_phi_local[dphi_map];
      double* phi_new_mapped = &ref_phi_new[phi_map];

      for (int g = 0; g < gss; g++)
        phi_new_mapped[g] += delta_phi_mapped[g];
    } // for dof
  }   // for cell
}
