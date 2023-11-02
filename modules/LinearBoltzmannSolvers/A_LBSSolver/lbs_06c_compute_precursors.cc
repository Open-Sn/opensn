#include "A_LBSSolver/lbs_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"


//###################################################################
/**Compute the steady state delayed neutron precursor concentrations.*/
void lbs::LBSSolver::ComputePrecursors()
{
  const size_t J = max_precursors_per_material_;

  precursor_new_local_.assign(precursor_new_local_.size(), 0.0);

  //================================================== Loop over cells
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& fe_values = unit_cell_matrices_[cell.local_id_];
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const double cell_volume = transport_view.Volume();

    //==================== Obtain xs
    const auto& xs = transport_view.XS();
    const auto& precursors = xs.Precursors();
    const auto& nu_delayed_sigma_f = xs.NuDelayedSigmaF();

    //======================================== Loop over precursors
    for (unsigned int j = 0; j < xs.NumPrecursors(); ++j)
    {
      size_t dof = cell.local_id_ * J + j;
      const auto& precursor = precursors[j];
      const double coeff = precursor.fractional_yield /
                           precursor.decay_constant;

      //=================================== Loop over nodes
      for (int i = 0; i < transport_view.NumNodes(); ++i)
      {
        const size_t uk_map = transport_view.MapDOF(i, 0, 0);
        const double node_V_fraction = fe_values.Vi_vectors[i]/cell_volume;

        //============================== Loop over groups
        for (unsigned int g = 0; g < groups_.size(); ++g)
          precursor_new_local_[dof] += coeff *
                                       nu_delayed_sigma_f[g] *
                                       phi_new_local_[uk_map + g] *
                                       node_V_fraction;
      }//for node i
    }//for precursor j

  }//for cell
}