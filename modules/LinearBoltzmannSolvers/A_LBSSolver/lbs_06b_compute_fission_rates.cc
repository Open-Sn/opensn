#include "A_LBSSolver/lbs_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
using namespace lbs;

//###################################################################
/**Compute the total fission production in the problem.
\author Zachary Hardy.*/
double
LBSSolver::ComputeFissionProduction(const std::vector<double>& phi)
{
  const int first_grp = groups_.front().id_;
  const int last_grp = groups_.back().id_;

  //============================================= Loop over local cells
  double local_production = 0.0;
  for (auto& cell : grid_ptr_->local_cells)
  {
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& cell_matrices = unit_cell_matrices_[cell.local_id_];

    //====================================== Obtain xs
    const auto& xs = transport_view.XS();
    const auto& F = xs.ProductionMatrix();
    const auto& nu_delayed_sigma_f = xs.NuDelayedSigmaF();

    if (not xs.IsFissionable()) continue;

    //====================================== Loop over nodes
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.Vi_vectors[i];

      //=============================== Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
      {
        const auto& prod = F[g];
        for (size_t gp = 0; gp <= last_grp; ++gp)
          local_production += prod[gp] * phi[uk_map + gp] * IntV_ShapeI;

        if (options_.use_precursors)
          for (unsigned int j = 0; j < xs.NumPrecursors(); ++j)
            local_production += nu_delayed_sigma_f[g] * phi[uk_map + g] * IntV_ShapeI;
      }
    } // for node
  }   // for cell

  //============================================= Allreduce global production
  double global_production = 0.0;
  MPI_Allreduce(&local_production,  // sendbuf
                &global_production, // recvbuf
                1,
                MPI_DOUBLE,     // count+datatype
                MPI_SUM,        // operation
                Chi::mpi.comm); // communicator

  return global_production;
}

//###################################################################
/**Computes the total fission rate in the problem.
\author Zachary Hardy.*/
double
LBSSolver::ComputeFissionRate(const std::vector<double>& phi)
{
  const int first_grp = groups_.front().id_;
  const int last_grp = groups_.back().id_;

  //============================================= Loop over local cells
  double local_fission_rate = 0.0;
  for (auto& cell : grid_ptr_->local_cells)
  {
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& cell_matrices = unit_cell_matrices_[cell.local_id_];

    //====================================== Obtain xs
    const auto& xs = transport_view.XS();
    const auto& sigma_f = xs.SigmaFission();

    // skip non-fissionable material
    if (not xs.IsFissionable()) continue;

    //====================================== Loop over nodes
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.Vi_vectors[i];

      //=============================== Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
        local_fission_rate += sigma_f[g] * phi[uk_map + g] * IntV_ShapeI;
    } // for node
  }   // for cell

  //============================================= Allreduce global production
  double global_fission_rate = 0.0;
  MPI_Allreduce(&local_fission_rate,  // sendbuf
                &global_fission_rate, // recvbuf
                1,
                MPI_DOUBLE,     // count+datatype
                MPI_SUM,        // operation
                Chi::mpi.comm); // communicator

  return global_fission_rate;
}