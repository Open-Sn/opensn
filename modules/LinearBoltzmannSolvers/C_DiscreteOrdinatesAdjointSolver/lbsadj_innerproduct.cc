#include "lbsadj_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

// ###################################################################
/**Computes the inner product of the flux and the material source.*/
double lbs::DiscreteOrdinatesAdjointSolver::ComputeInnerProduct()
{
  double local_integral = 0.0;

  //============================================= Material sources
  for (const auto& cell : grid_ptr_->local_cells)
  {
    if (matid_to_src_map_.count(cell.material_id_) == 0) continue; //Skip if no src

    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& source = matid_to_src_map_[cell.material_id_];
    const auto& fe_values = unit_cell_matrices_[cell.local_id_];

    for (const auto& group : groups_)
    {
      const int g = group.id_;
      const double Q = source->source_value_g_[g];

      if (Q > 0.0)
      {
        const int num_nodes = transport_view.NumNodes();
        for (int i = 0; i < num_nodes; ++i)
        {
          const size_t dof_map = transport_view.MapDOF(i, 0, g); //unknown map

          const double phi = phi_old_local_[dof_map];

          local_integral += Q * phi * fe_values.Vi_vectors[i];
        }//for node
      }//check source value >0
    }//for group
  }//for cell

  //============================================= Point sources
  for (const auto& point_source : point_sources_)
  {
    const auto& info_list = point_source.ContainingCellsInfo();
    for (const auto& info : info_list)
    {
      const auto& cell = grid_ptr_->local_cells[info.cell_local_id];
      const auto& transport_view = cell_transport_views_[cell.local_id_];
      const auto& source_strength = point_source.Strength();
      const auto& shape_values = info.shape_values;

      for (const auto& group : groups_)
      {
        const int g = group.id_;
        const double S = source_strength[g] * info.volume_weight;

        if (S > 0.0)
        {
          const int num_nodes = transport_view.NumNodes();
          for (int i = 0; i < num_nodes; ++i)
          {
            const size_t dof_map = transport_view.MapDOF(i, 0, g); //unknown map

            const double phi_i = phi_old_local_[dof_map];

            local_integral += S * phi_i * shape_values[i];
          }//for node
        }//check source value >0
      }//for group
    }//for cell
  }//for point source

  double global_integral = 0.0;

  MPI_Allreduce(&local_integral,     //sendbuf
                &global_integral,    //recvbuf
                1, MPI_DOUBLE,       //count, datatype
                MPI_SUM,             //op
                Chi::mpi.comm);     //comm

  return global_integral;
}