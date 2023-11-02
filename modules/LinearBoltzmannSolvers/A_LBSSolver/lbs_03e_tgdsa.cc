#include "A_LBSSolver/lbs_solver.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "A_LBSSolver/Acceleration/diffusion_mip.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

// ###################################################################
/**Initializes the Within-Group DSA solver. */
void lbs::LBSSolver::InitTGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa_)
  {
    //=========================================== Make UnknownManager
    const auto& uk_man = discretization_->UNITARY_UNKNOWN_MANAGER;

    //=========================================== Make boundary conditions
    auto bcs = acceleration::TranslateBCs(sweep_boundaries_);

    //=========================================== Make TwoGridInfo
     for (const auto& mat_id_xs_pair : matid_to_xs_map_)
    {
       const auto& mat_id = mat_id_xs_pair.first;
       const auto& xs     = mat_id_xs_pair.second;

       acceleration::TwoGridCollapsedInfo tginfo =
         MakeTwoGridCollapsedInfo(*xs,
         acceleration::EnergyCollapseScheme::JFULL);

       groupset.tg_acceleration_info_.map_mat_id_2_tginfo.insert(
         std::make_pair(mat_id, std::move(tginfo)));
     }

    //=========================================== Make xs map
     typedef lbs::acceleration::Multigroup_D_and_sigR MGXS;
     typedef std::map<int, MGXS> MatID2MGDXSMap;
     MatID2MGDXSMap matid_2_mgxs_map;
     for (const auto& matid_xs_pair : matid_to_xs_map_)
    {
       const auto& mat_id = matid_xs_pair.first;

       const auto& tg_info =
         groupset.tg_acceleration_info_.map_mat_id_2_tginfo.at(mat_id);

       matid_2_mgxs_map.insert(
         std::make_pair(mat_id, MGXS{{tg_info.collapsed_D},
                                     {tg_info.collapsed_sig_a}}));
     }

    //=========================================== Create solver
    const auto& sdm = *discretization_;

    auto solver = std::make_shared<acceleration::DiffusionMIPSolver>(
      std::string(TextName() + "_TGDSA"),
      sdm,
      uk_man,
      bcs,
      matid_2_mgxs_map,
      unit_cell_matrices_,
      true); // verbosity

    solver->options.residual_tolerance = groupset.tgdsa_tol_;
    solver->options.max_iters = groupset.tgdsa_max_iters_;
    solver->options.verbose = groupset.tgdsa_verbose_;
    solver->options.additional_options_string = groupset.tgdsa_string_;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);

    solver->AssembleAand_b(dummy_rhs);

    groupset.tgdsa_solver_ = solver;
  }
}

// ###################################################################
/**Cleans up memory consuming items. */
void lbs::LBSSolver::CleanUpTGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa_) groupset.tgdsa_solver_ = nullptr;
}

// ###################################################################
/**Assembles a delta-phi vector on the first moment.*/
void lbs::LBSSolver::AssembleTGDSADeltaPhiVector(
  const LBSGroupset& groupset,
  const std::vector<double>& phi_in,
  std::vector<double>& delta_phi_local)
{
  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  delta_phi_local.clear();
  delta_phi_local.assign(local_node_count_, 0.0);

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& S = matid_to_xs_map_[cell.material_id_]->TransferMatrix(0);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

      double& delta_phi_mapped = delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g = 0; g < gss; ++g)
      {
        double R_g = 0.0;
        for (const auto& [row_g, gprime, sigma_sm] : S.Row(gsi + g))
          if (gprime >= gsi and gprime != (gsi + g))
            R_g += sigma_sm * phi_in_mapped[gprime];

        delta_phi_mapped += R_g;
      } // for g
    }   // for node
  }     // for cell
}

// ###################################################################
/**DAssembles a delta-phi vector on the first moment.*/
void lbs::LBSSolver::DisAssembleTGDSADeltaPhiVector(
  const LBSGroupset& groupset,
  const std::vector<double>& delta_phi_local,
  std::vector<double>& ref_phi_new)
{
  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  const auto& map_mat_id_2_tginfo =
    groupset.tg_acceleration_info_.map_mat_id_2_tginfo;

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto& xi_g = map_mat_id_2_tginfo.at(cell.material_id_).spectrum;

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double delta_phi_mapped = delta_phi_local[dphi_map];
      double* phi_new_mapped = &ref_phi_new[phi_map];

      for (int g = 0; g < gss; ++g)
        phi_new_mapped[g] += delta_phi_mapped * xi_g[gsi + g];
    } // for dof
  }   // for cell
}