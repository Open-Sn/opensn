#include "lbs_discrete_ordinates_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include <iomanip>

//###################################################################
/**Zeroes all the outflow data-structures required to compute
 * balance.*/
void
lbs::DiscreteOrdinatesSolver::ZeroOutflowBalanceVars(LBSGroupset& groupset)
{
  for (auto& cell_transport_view : cell_transport_views_)
    for (auto& group : groupset.groups_)
      cell_transport_view.ZeroOutflow(group.id_);
}

//###################################################################
/**Compute balance.*/
void
lbs::DiscreteOrdinatesSolver::ComputeBalance()
{
  Chi::mpi.Barrier();
  Chi::log.Log() << "\n********** Computing balance\n";

  //======================================== Get material source
  // This is done using the SetSource routine
  // because it allows a lot of flexibility.
  auto mat_src = phi_old_local_;
  mat_src.assign(mat_src.size(), 0.0);
  for (auto& groupset : groupsets_)
  {
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    active_set_source_function_(groupset,
                                q_moments_local_,
                                PhiOldLocal(),
                                APPLY_FIXED_SOURCES | APPLY_AGS_FISSION_SOURCES |
                                  APPLY_WGS_FISSION_SOURCES);
    LBSSolver::GSScopedCopyPrimarySTLvectors(groupset, q_moments_local_, mat_src);
  }

  //======================================== Compute absorption, material-source
  //                                         and in-flow
  double local_out_flow = 0.0;
  double local_in_flow = 0.0;
  double local_absorption = 0.0;
  double local_production = 0.0;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = discretization_->GetCellMapping(cell);
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& fe_intgrl_values = unit_cell_matrices_[cell.local_id_];
    const size_t num_nodes = transport_view.NumNodes();
    const auto& IntV_shapeI = fe_intgrl_values.Vi_vectors;
    const auto& IntS_shapeI = fe_intgrl_values.face_Si_vectors;

    //====================================== Inflow
    // This is essentially an integration over
    // all faces, all angles, and all groups.
    // Only the cosines that are negative are
    // added to the integral.
    for (int f = 0; f < cell.faces_.size(); ++f)
    {
      const auto& face = cell.faces_[f];

      for (const auto& groupset : groupsets_)
      {
        for (int n = 0; n < groupset.quadrature_->omegas_.size(); ++n)
        {
          const auto& omega = groupset.quadrature_->omegas_[n];
          const double wt = groupset.quadrature_->weights_[n];
          const double mu = omega.Dot(face.normal_);

          if (mu < 0.0 and (not face.has_neighbor_)) // mu<0 and bndry
          {
            const auto& bndry = sweep_boundaries_[face.neighbor_id_];
            for (int fi = 0; fi < face.vertex_ids_.size(); ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const auto& IntFi_shapeI = IntS_shapeI[f][i];

              for (const auto& group : groupset.groups_)
              {
                const int g = group.id_;
                const double psi = *bndry->HeterogeneousPsiIncoming(cell.local_id_, f, fi, n, g, 0);
                local_in_flow -= mu * wt * psi * IntFi_shapeI;
              } // for g
            }   // for fi
          }     // if bndry
        }       // for n
      }         // for groupset
    }           // for f

    //====================================== Outflow
    // The group-wise outflow was determined
    // during a solve so here we just
    // consolidate it.
    for (int g = 0; g < num_groups_; ++g)
      local_out_flow += transport_view.GetOutflow(g);

    //====================================== Absorption and Src
    // Isotropic flux based absorption and source
    const auto& xs = transport_view.XS();
    const auto& sigma_a = xs.SigmaAbsorption();
    for (int i = 0; i < num_nodes; ++i)
      for (int g = 0; g < num_groups_; ++g)
      {
        size_t imap = transport_view.MapDOF(i, 0, g);
        double phi_0g = phi_old_local_[imap];
        double q_0g = mat_src[imap];

        local_absorption += sigma_a[g] * phi_0g * IntV_shapeI[i];
        local_production += q_0g * IntV_shapeI[i];
      } // for g
  }     // for cell

  //======================================== Consolidate local balances
  double local_balance = local_production + local_in_flow - local_absorption - local_out_flow;
  double local_gain = local_production + local_in_flow;

  std::vector<double> local_balance_table = {
    local_absorption, local_production, local_in_flow, local_out_flow, local_balance, local_gain};
  size_t table_size = local_balance_table.size();

  std::vector<double> globl_balance_table(table_size, 0.0);

  MPI_Allreduce(local_balance_table.data(), // sendbuf
                globl_balance_table.data(), // recvbuf
                table_size,
                MPI_DOUBLE,     // count + datatype
                MPI_SUM,        // operation
                Chi::mpi.comm); // communicator

  double globl_absorption = globl_balance_table.at(0);
  double globl_production = globl_balance_table.at(1);
  double globl_in_flow = globl_balance_table.at(2);
  double globl_out_flow = globl_balance_table.at(3);
  double globl_balance = globl_balance_table.at(4);
  double globl_gain = globl_balance_table.at(5);

  Chi::log.Log() << "Balance table:\n"
                 << std::setprecision(5) << std::scientific
                 << " Absorption rate              = " << globl_absorption << "\n"
                 << " Production rate              = " << globl_production << "\n"
                 << " In-flow rate                 = " << globl_in_flow << "\n"
                 << " Out-flow rate                = " << globl_out_flow << "\n"
                 << " Net Gain (In-flow + sources) = " << globl_gain << "\n"
                 << " Net Balance                  = " << globl_balance << "\n"
                 << " (Net Balance)/(Net Gain)     = " << globl_balance / globl_gain << "\n";

  Chi::log.Log() << "\n********** Done computing balance\n";

  Chi::mpi.Barrier();
}