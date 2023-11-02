#include "lbs_discrete_ordinates_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

/**Computes the angular flux based leakage from boundary surfaces.
\param groupset_id The groupset for which to compute the leakage.
\param boundary_id uint64_t The boundary-id for which to perform the
integration.

\return The leakage as a value.
*/
std::vector<double>
lbs::DiscreteOrdinatesSolver::ComputeLeakage(const int groupset_id,
                                             const uint64_t boundary_id) const
{
  const std::string fname = "lbs::SteadySolver::ComputeLeakage";

  //================================================== Perform checks
  if (groupset_id < 0 or groupset_id >= groupsets_.size())
    throw std::invalid_argument(fname + ": Invalid groupset_id specified.");

  if (not options_.save_angular_flux)
    throw std::logic_error(fname + ": Requires options.save_angular_flux to be"
                                   " true.");

  //================================================== Get info
  const auto& sdm = *discretization_;
  const auto& groupset = groupsets_.at(groupset_id);
  const auto& psi_uk_man = groupset.psi_uk_man_;
  const auto& quadrature = groupset.quadrature_;

  const size_t num_angles = quadrature->omegas_.size();

  const int gsi = groupset.groups_.front().id_;
  const int gsf = groupset.groups_.back().id_;
  const int gs_num_groups = gsf + 1 - gsi;

  //================================================== Start integration
  std::vector<double> local_leakage(gs_num_groups, 0.0);
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& fe_values = unit_cell_matrices_[cell.local_id_];

    size_t f = 0;
    for (const auto& face : cell.faces_)
    {
      if (not face.has_neighbor_ and face.neighbor_id_ == boundary_id)
      {
        const auto& IntF_shapeI = fe_values.face_Si_vectors[f];
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          for (size_t n = 0; n < num_angles; ++n)
          {
            const auto& omega = quadrature->omegas_[n];
            const auto& weight = quadrature->weights_[n];
            const double mu = omega.Dot(face.normal_);
            if (mu > 0.0)
            {
              for (int gi = 0; gi < gs_num_groups; ++gi)
              {
                const int g = gi + gsi;
                const int64_t imap = sdm.MapDOFLocal(cell, i, psi_uk_man, n, g);

                const double psi = psi_new_local_[groupset_id][imap];

                local_leakage[gi] += weight * mu * psi * IntF_shapeI[i];
              } // for g
            }   // outgoing
          }     // for n
        }       // for face node
      }         // if right bndry
      ++f;
    } // for face
  }   // for cell

  std::vector<double> global_leakage(gs_num_groups, 0.0);
  MPI_Allreduce(local_leakage.data(),  // sendbuf
                global_leakage.data(), // recvbuf,
                gs_num_groups,
                MPI_DOUBLE,     // count+datatype
                MPI_SUM,        // operation
                Chi::mpi.comm); // comm

  return global_leakage;
}