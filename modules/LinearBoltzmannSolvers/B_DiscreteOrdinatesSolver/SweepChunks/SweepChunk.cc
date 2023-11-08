#include "modules/LinearBoltzmannSolvers/B_DiscreteOrdinatesSolver/SweepChunks/SweepChunk.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"
#include "framework/math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/logging/chi_log_exceptions.h"

#define scint static_cast<int>

namespace lbs
{

SweepChunk::SweepChunk(std::vector<double>& destination_phi,
                       std::vector<double>& destination_psi,
                       const chi_mesh::MeshContinuum& grid,
                       const chi_math::SpatialDiscretization& discretization,
                       const std::vector<UnitCellMatrices>& unit_cell_matrices,
                       std::vector<lbs::CellLBSView>& cell_transport_views,
                       const std::vector<double>& source_moments,
                       const LBSGroupset& groupset,
                       const std::map<int, XSPtr>& xs,
                       int num_moments,
                       int max_num_cell_dofs,
                       std::unique_ptr<SweepDependencyInterface> sweep_dependency_interface_ptr)
  : chi_mesh::sweep_management::SweepChunk(destination_phi, destination_psi),
    grid_(grid),
    grid_fe_view_(discretization),
    unit_cell_matrices_(unit_cell_matrices),
    grid_transport_view_(cell_transport_views),
    q_moments_(source_moments),
    groupset_(groupset),
    xs_(xs),
    num_moments_(num_moments),
    save_angular_flux_(!destination_psi.empty()),
    sweep_dependency_interface_ptr_(std::move(sweep_dependency_interface_ptr)),
    sweep_dependency_interface_(*sweep_dependency_interface_ptr_),
    groupset_angle_group_stride_(groupset_.psi_uk_man_.NumberOfUnknowns() *
                                 groupset_.groups_.size()),
    groupset_group_stride_(groupset_.groups_.size())
{
  Amat_.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
  Atemp_.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
  b_.resize(groupset.groups_.size(), std::vector<double>(max_num_cell_dofs, 0.0));
  source_.resize(max_num_cell_dofs, 0.0);

  sweep_dependency_interface_.groupset_angle_group_stride_ = groupset_angle_group_stride_;
  sweep_dependency_interface_.groupset_group_stride_ = groupset_group_stride_;
}

// ##################################################################
/**Registers a kernel as a named callback function*/
void
SweepChunk::RegisterKernel(const std::string& name, CallbackFunction function)
{
  ChiInvalidArgumentIf(kernels_.count(name) > 0,
                       "Attempting to register kernel with name \"" + name +
                         "\" but the kernel already exists.");

  kernels_[name] = std::move(function);
}

// ##################################################################
/**Returns a kernel if the given name exists.*/
SweepChunk::CallbackFunction
SweepChunk::Kernel(const std::string& name) const
{
  ChiInvalidArgumentIf(kernels_.count(name) == 0,
                       "No register kernel with name \"" + name + "\" found");
  return kernels_.at(name);
}

// ##################################################################
/**Executes the supplied kernels list.*/
void
SweepChunk::ExecuteKernels(const std::vector<CallbackFunction>& kernels)
{
  for (auto& kernel : kernels)
    kernel();
}

// ##################################################################
/**Operations when outgoing fluxes are handled including passing
 * face angular fluxes downstream and computing
 * balance parameters (i.e. outflow)
 * */
void
SweepChunk::OutgoingSurfaceOperations()
{
  const size_t f = sweep_dependency_interface_.current_face_idx_;
  const auto& IntF_shapeI = (*IntS_shapeI_)[f];
  const double mu = face_mu_values_[f];
  const double wt = direction_qweight_;

  const bool on_boundary = sweep_dependency_interface_.on_boundary_;
  const bool is_reflecting_boundary = sweep_dependency_interface_.is_reflecting_bndry_;

  const size_t num_face_nodes = cell_mapping_->NumFaceNodes(f);
  for (int fi = 0; fi < num_face_nodes; ++fi)
  {
    const int i = cell_mapping_->MapFaceNode(f, fi);

    double* psi = sweep_dependency_interface_.GetDownwindPsi(fi);

    if (psi != nullptr)
      if (not on_boundary or is_reflecting_boundary)
        for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
          psi[gsg] = b_[gsg][i];
    if (on_boundary and not is_reflecting_boundary)
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        cell_transport_view_->AddOutflow(gs_gi_ + gsg, wt * mu * b_[gsg][i] * IntF_shapeI[i]);

  } // for fi
}

// ##################################################################
/**Assembles the volumetric gradient term.*/
void
SweepChunk::KernelFEMVolumetricGradientTerm()
{
  const auto& G = *G_;

  for (int i = 0; i < cell_num_nodes_; ++i)
    for (int j = 0; j < cell_num_nodes_; ++j)
      Amat_[i][j] = omega_.Dot(G[i][j]);
}

// ##################################################################
/**Performs the integral over the surface of a face.*/
void
SweepChunk::KernelFEMUpwindSurfaceIntegrals()
{
  const size_t f = sweep_dependency_interface_.current_face_idx_;
  const auto& M_surf_f = (*M_surf_)[f];
  const double mu = face_mu_values_[f];
  const size_t num_face_nodes = sweep_dependency_interface_.num_face_nodes_;
  for (int fi = 0; fi < num_face_nodes; ++fi)
  {
    const int i = cell_mapping_->MapFaceNode(f, fi);
    for (int fj = 0; fj < num_face_nodes; ++fj)
    {
      const int j = cell_mapping_->MapFaceNode(f, fj);

      const double* psi = sweep_dependency_interface_.GetUpwindPsi(fj);

      const double mu_Nij = -mu * M_surf_f[i][j];
      Amat_[i][j] += mu_Nij;

      if (psi == nullptr) continue;

      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        b_[gsg][i] += psi[gsg] * mu_Nij;
    } // for face node j
  }   // for face node i
}

// ##################################################################
/**Assembles angular sources and applies the mass matrix terms.*/
void
SweepChunk::KernelFEMSTDMassTerms()
{
  const auto& M = *M_;
  const auto& m2d_op = groupset_.quadrature_->GetMomentToDiscreteOperator();

  // ============================= Contribute source moments
  // q = M_n^T * q_moms
  for (int i = 0; i < cell_num_nodes_; ++i)
  {
    double temp_src = 0.0;
    for (int m = 0; m < num_moments_; ++m)
    {
      const size_t ir = cell_transport_view_->MapDOF(i, m, scint(g_));
      temp_src += m2d_op[m][direction_num_] * q_moments_[ir];
    } // for m
    source_[i] = temp_src;
  } // for i

  // ============================= Mass Matrix and Source
  // Atemp  = Amat + sigma_tgr * M
  // b     += M * q
  for (int i = 0; i < cell_num_nodes_; ++i)
  {
    double temp = 0.0;
    for (int j = 0; j < cell_num_nodes_; ++j)
    {
      const double Mij = M[i][j];
      Atemp_[i][j] = Amat_[i][j] + Mij * sigma_tg_;
      temp += Mij * source_[j];
    } // for j
    b_[gsg_][i] += temp;
  } // for i
}

// ##################################################################
/**Adds a single direction's contribution to the moment integrals.*/
void
SweepChunk::KernelPhiUpdate()
{
  const auto& d2m_op = groupset_.quadrature_->GetDiscreteToMomentOperator();

  auto& output_phi = GetDestinationPhi();

  for (int m = 0; m < num_moments_; ++m)
  {
    const double wn_d2m = d2m_op[m][direction_num_];
    for (int i = 0; i < cell_num_nodes_; ++i)
    {
      const size_t ir = cell_transport_view_->MapDOF(i, m, gs_gi_);
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        output_phi[ir + gsg] += wn_d2m * b_[gsg][i];
    }
  }
}

// ##################################################################
/**Updates angular fluxes.*/
void
SweepChunk::KernelPsiUpdate()
{
  if (not save_angular_flux_) return;

  auto& output_psi = GetDestinationPsi();
  double* cell_psi_data =
    &output_psi[grid_fe_view_.MapDOFLocal(*cell_, 0, groupset_.psi_uk_man_, 0, 0)];

  for (size_t i = 0; i < cell_num_nodes_; ++i)
  {
    const size_t imap =
      i * groupset_angle_group_stride_ + direction_num_ * groupset_group_stride_ + gs_ss_begin_;
    for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
      cell_psi_data[imap + gsg] = b_[gsg][i];
  } // for i
}

// ##################################################################
/**Sets data for the current incoming face.*/
void
SweepDependencyInterface::SetupIncomingFace(
  int face_id, size_t num_face_nodes, uint64_t neighbor_id, bool on_local_face, bool on_boundary)
{
  current_face_idx_ = face_id;
  num_face_nodes_ = num_face_nodes;
  neighbor_id_ = neighbor_id;
  on_local_face_ = on_local_face;
  on_boundary_ = on_boundary;
}

// ##################################################################
/**Sets data for the current outgoing face.*/
void
SweepDependencyInterface::SetupOutgoingFace(int face_id,
                                            size_t num_face_nodes,
                                            uint64_t neighbor_id,
                                            bool on_local_face,
                                            bool on_boundary,
                                            int locality)
{
  current_face_idx_ = face_id;
  num_face_nodes_ = num_face_nodes;
  neighbor_id_ = neighbor_id;
  face_locality_ = locality;
  on_local_face_ = on_local_face;
  on_boundary_ = on_boundary;

  is_reflecting_bndry_ =
    (on_boundary_ and angle_set_->GetBoundaries()[neighbor_id_]->IsReflecting());
}

} // namespace lbs
