#include "lbsadj_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

/**Main execution function.*/
void lbs::DiscreteOrdinatesAdjointSolver::Execute()
{
  const std::string fname = __FUNCTION__;

  primary_ags_solver_->Setup();
  primary_ags_solver_->Solve();

  //============================================= Apply post processing
  Chi::log.Log() << "LBAdjointSolver: post-processing.";
  std::set<int> set_group_numbers;
  for (const auto& groupset : groupsets_)
    for (const auto& group : groupset.groups_)
      set_group_numbers.insert(group.id_);

  const auto& m_to_ell_em_map =
    groupsets_.front().quadrature_->GetMomentToHarmonicsIndexMap();

  //============================================= Reorient phi-moments for reverse
  //                                              angle
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_view = cell_transport_views_[cell.local_id_];
    const int num_nodes = cell_view.NumNodes();

    for (int i = 0; i < num_nodes; ++i)
    {
      for (int m = 0; m < num_moments_; ++m)
      {
        const auto& ell = m_to_ell_em_map[m].ell;

        size_t dof_map_g0 = cell_view.MapDOF(i, m, 0); //unknown map

        for (int g : set_group_numbers)
          phi_old_local_[dof_map_g0 + g] *= pow(-1.0, ell);
      }//for moment
    }//node i
  }//for cell

  UpdateFieldFunctions();

}