#include "lbs_solver.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

namespace lbs
{

// ###################################################################
/**Copy relevant section of phi_old to the field functions.*/
void
LBSSolver::UpdateFieldFunctions()
{
  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  //======================================== Update flux moments
  for (const auto& [g_and_m, ff_index] : phi_field_functions_local_map_)
  {
    const size_t g = g_and_m.first;
    const size_t m = g_and_m.second;

    std::vector<double> data_vector_local(local_node_count_, 0.0);

    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imapA = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
        const int64_t imapB = sdm.MapDOFLocal(cell, i);

        data_vector_local[imapB] = phi_old_local_[imapA];
      } // for node
    }   // for cell

    auto& ff_ptr = field_functions_.at(ff_index);
    ff_ptr->UpdateFieldVector(data_vector_local);
  }
  // for (size_t g = 0; g < groups_.size(); ++g)
  //{
  //   for (size_t m = 0; m < num_moments_; ++m)
  //   {
  //     std::vector<double> data_vector_local(local_node_count_, 0.0);
  //
  //     for (const auto& cell : grid_ptr_->local_cells)
  //     {
  //       const auto& cell_mapping = sdm.GetCellMapping(cell);
  //       const size_t num_nodes = cell_mapping.NumNodes();
  //
  //       for (size_t i = 0; i < num_nodes; ++i)
  //       {
  //         const int64_t imapA = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
  //         const int64_t imapB = sdm.MapDOFLocal(cell, i);
  //
  //         data_vector_local[imapB] = phi_old_local_[imapA];
  //       } // for node
  //     }   // for cell
  //
  //     ChiLogicalErrorIf(phi_field_functions_local_map_.count({g, m}) == 0,
  //                       "Update error for phi based field functions");
  //
  //     const size_t ff_index = phi_field_functions_local_map_.at({g, m});
  //
  //     auto& ff_ptr = field_functions_.at(ff_index);
  //     ff_ptr->UpdateFieldVector(data_vector_local);
  //   } // for m
  // }   // for g

  //======================================== Update power generation
  //                                         if enabled
  if (options_.power_field_function_on)
  {
    std::vector<double> data_vector_local(local_node_count_, 0.0);

    double local_total_power = 0.0;
    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      const auto& Vi = unit_cell_matrices_[cell.local_id_].Vi_vectors;

      const auto& xs = matid_to_xs_map_.at(cell.material_id_);

      if (not xs->IsFissionable()) continue;

      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imapA = sdm.MapDOFLocal(cell, i);
        const int64_t imapB = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

        double nodal_power = 0.0;
        for (size_t g = 0; g < groups_.size(); ++g)
        {
          const double sigma_fg = xs->SigmaFission()[g];
          // const double kappa_g = xs->Kappa()[g];
          const double kappa_g = options_.power_default_kappa;

          nodal_power += kappa_g * sigma_fg * phi_old_local_[imapB + g];
        } // for g

        data_vector_local[imapA] = nodal_power;
        local_total_power += nodal_power * Vi[i];
      } // for node
    }   // for cell

    if (options_.power_normalization > 0.0)
    {
      double globl_total_power;
      MPI_Allreduce(&local_total_power, // sendbuf
                    &globl_total_power, // recvbuf
                    1,
                    MPI_DOUBLE,     // count + datatype
                    MPI_SUM,        // operation
                    Chi::mpi.comm); // communicator

      chi_math::Scale(data_vector_local, options_.power_normalization / globl_total_power);
    }

    const size_t ff_index = power_gen_fieldfunc_local_handle_;

    auto& ff_ptr = field_functions_.at(ff_index);
    ff_ptr->UpdateFieldVector(data_vector_local);

  } // if power enabled
}

// ###################################################################
/**Sets the internal phi vector to the value in the associated
field function.*/
void
LBSSolver::SetPhiFromFieldFunctions(PhiSTLOption which_phi,
                                    const std::vector<size_t>& m_indices,
                                    const std::vector<size_t>& g_indices)
{
  std::vector<size_t> m_ids_to_copy = m_indices;
  std::vector<size_t> g_ids_to_copy = g_indices;
  if (m_indices.empty())
    for (size_t m = 0; m < num_moments_; ++m)
      m_ids_to_copy.push_back(m);
  if (g_ids_to_copy.empty())
    for (size_t g = 0; g < num_groups_; ++g)
      g_ids_to_copy.push_back(g);

  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  for (const size_t m : m_ids_to_copy)
  {
    for (const size_t g : g_ids_to_copy)
    {
      const size_t ff_index = phi_field_functions_local_map_.at({g, m});
      auto& ff_ptr = field_functions_.at(ff_index);
      auto& ff_data = ff_ptr->FieldVector();

      for (const auto& cell : grid_ptr_->local_cells)
      {
        const auto& cell_mapping = sdm.GetCellMapping(cell);
        const size_t num_nodes = cell_mapping.NumNodes();

        for (size_t i = 0; i < num_nodes; ++i)
        {
          const int64_t imapA = sdm.MapDOFLocal(cell, i);
          const int64_t imapB = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);

          if (which_phi == PhiSTLOption::PHI_OLD) phi_old_local_[imapB] = ff_data[imapA];
          else if (which_phi == PhiSTLOption::PHI_NEW)
            phi_new_local_[imapB] = ff_data[imapA];
        } // for node
      }   // for cell
    }     // for g
  }       // for m
}

} // namespace lbs
