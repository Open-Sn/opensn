#include "modules/LinearBoltzmannSolvers/B_DiscreteOrdinatesSolver/SweepChunks/AAH_SweepChunk.h"

#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "framework/mesh/SweepUtilities/FLUDS/AAH_FLUDS.h"

#define scint static_cast<int>

namespace lbs
{

AAH_SweepChunk::AAH_SweepChunk(const chi_mesh::MeshContinuum& grid,
                               const chi_math::SpatialDiscretization& discretization,
                               const std::vector<UnitCellMatrices>& unit_cell_matrices,
                               std::vector<lbs::CellLBSView>& cell_transport_views,
                               std::vector<double>& destination_phi,
                               std::vector<double>& destination_psi,
                               const std::vector<double>& source_moments,
                               const LBSGroupset& groupset,
                               const std::map<int, XSPtr>& xs,
                               int num_moments,
                               int max_num_cell_dofs)
  : SweepChunk(destination_phi,
               destination_psi,
               grid,
               discretization,
               unit_cell_matrices,
               cell_transport_views,
               source_moments,
               groupset,
               xs,
               num_moments,
               max_num_cell_dofs,
               std::make_unique<AAH_SweepDependencyInterface>())
{
  // Register kernels
  RegisterKernel("FEMVolumetricGradTerm",
                 std::bind(&SweepChunk::KernelFEMVolumetricGradientTerm, this));
  RegisterKernel("FEMUpwindSurfaceIntegrals",
                 std::bind(&SweepChunk::KernelFEMUpwindSurfaceIntegrals, this));
  RegisterKernel("FEMSSTDMassTerms", std::bind(&SweepChunk::KernelFEMSTDMassTerms, this));
  RegisterKernel("KernelPhiUpdate", std::bind(&SweepChunk::KernelPhiUpdate, this));
  RegisterKernel("KernelPsiUpdate", std::bind(&SweepChunk::KernelPsiUpdate, this));

  // Setup callbacks
  cell_data_callbacks_ = {};

  direction_data_callbacks_and_kernels_ = {Kernel("FEMVolumetricGradTerm")};

  surface_integral_kernels_ = {Kernel("FEMUpwindSurfaceIntegrals")};

  mass_term_kernels_ = {Kernel("FEMSSTDMassTerms")};

  flux_update_kernels_ = {Kernel("KernelPhiUpdate"), Kernel("KernelPsiUpdate")};

  post_cell_dir_sweep_callbacks_ = {};
}

void
AAH_SweepChunk::Sweep(chi_mesh::sweep_management::AngleSet& angle_set)
{
  const chi::SubSetInfo& grp_ss_info = groupset_.grp_subset_infos_[angle_set.GetRefGroupSubset()];

  gs_ss_size_ = grp_ss_info.ss_size;
  gs_ss_begin_ = grp_ss_info.ss_begin;
  gs_gi_ = groupset_.groups_[gs_ss_begin_].id_;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  sweep_dependency_interface_.angle_set_ = &angle_set;
  sweep_dependency_interface_.surface_source_active_ = IsSurfaceSourceActive();
  sweep_dependency_interface_.gs_ss_begin_ = gs_ss_begin_;
  sweep_dependency_interface_.gs_gi_ = gs_gi_;

  auto& aah_sweep_depinterf =
    dynamic_cast<AAH_SweepDependencyInterface&>(sweep_dependency_interface_);
  aah_sweep_depinterf.fluds_ =
    &dynamic_cast<chi_mesh::sweep_management::AAH_FLUDS&>(angle_set.GetFLUDS());

  // Loop over each cell
  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetSPLS().item_id;
  const size_t num_spls = spls.size();
  for (size_t spls_index = 0; spls_index < num_spls; ++spls_index)
  {
    cell_local_id_ = spls[spls_index];
    cell_ = &grid_.local_cells[cell_local_id_];
    sweep_dependency_interface_.cell_ptr_ = cell_;
    sweep_dependency_interface_.cell_local_id_ = cell_local_id_;
    cell_mapping_ = &grid_fe_view_.GetCellMapping(*cell_);
    cell_transport_view_ = &grid_transport_view_[cell_->local_id_];

    using namespace chi_mesh::sweep_management;
    const auto& face_orientations = spds.CellFaceOrientations()[cell_local_id_];

    cell_num_faces_ = cell_->faces_.size();
    cell_num_nodes_ = cell_mapping_->NumNodes();
    const auto& sigma_t = xs_.at(cell_->material_id_)->SigmaTotal();

    aah_sweep_depinterf.spls_index = spls_index;

    // Get Cell matrices
    const auto& fe_intgrl_values = unit_cell_matrices_[cell_local_id_];
    G_ = &fe_intgrl_values.G_matrix;
    M_ = &fe_intgrl_values.M_matrix;
    M_surf_ = &fe_intgrl_values.face_M_matrices;
    IntS_shapeI_ = &fe_intgrl_values.face_Si_vectors;

    for (auto& callback : cell_data_callbacks_)
      callback();

    // Loop over angles in set
    const int ni_deploc_face_counter = deploc_face_counter;
    const int ni_preloc_face_counter = preloc_face_counter;

    // as = angle set
    // ss = subset
    const std::vector<size_t>& as_angle_indices = angle_set.GetAngleIndices();
    const size_t as_num_angles = as_angle_indices.size();
    for (size_t as_ss_idx = 0; as_ss_idx < as_num_angles; ++as_ss_idx)
    {
      direction_num_ = as_angle_indices[as_ss_idx];
      omega_ = groupset_.quadrature_->omegas_[direction_num_];
      direction_qweight_ = groupset_.quadrature_->weights_[direction_num_];

      sweep_dependency_interface_.angle_set_index_ = as_ss_idx;
      sweep_dependency_interface_.angle_num_ = direction_num_;

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      // Reset right-handside
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        b_[gsg].assign(cell_num_nodes_, 0.0);

      ExecuteKernels(direction_data_callbacks_and_kernels_);

      // Upwinding structure
      aah_sweep_depinterf.in_face_counter = 0;
      aah_sweep_depinterf.preloc_face_counter = 0;
      aah_sweep_depinterf.out_face_counter = 0;
      aah_sweep_depinterf.deploc_face_counter = 0;

      // Update face orientations
      face_mu_values_.assign(cell_num_faces_, 0.0);
      for (int f = 0; f < cell_num_faces_; ++f)
        face_mu_values_[f] = omega_.Dot(cell_->faces_[f].normal_);

      // Surface integrals
      int in_face_counter = -1;
      for (int f = 0; f < cell_num_faces_; ++f)
      {
        const auto& face = cell_->faces_[f];

        if (face_orientations[f] != FaceOrientation::INCOMING) continue;

        const bool local = cell_transport_view_->IsFaceLocal(f);
        const bool boundary = not face.has_neighbor_;

        if (local) ++in_face_counter;
        else if (not boundary)
          ++preloc_face_counter;

        sweep_dependency_interface_.SetupIncomingFace(
          f, cell_mapping_->NumFaceNodes(f), face.neighbor_id_, local, boundary);

        aah_sweep_depinterf.in_face_counter = in_face_counter;
        aah_sweep_depinterf.preloc_face_counter = preloc_face_counter;

        // IntSf_mu_psi_Mij_dA
        ExecuteKernels(surface_integral_kernels_);
      } // for f

      // Looping over groups, assembling mass terms
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
      {
        g_ = gs_gi_ + gsg;
        gsg_ = gsg;
        sigma_tg_ = sigma_t[g_];

        ExecuteKernels(mass_term_kernels_);

        // Solve system
        chi_math::GaussElimination(Atemp_, b_[gsg], scint(cell_num_nodes_));
      }

      // Flux updates
      ExecuteKernels(flux_update_kernels_);

      // Perform outgoing surface operations
      int out_face_counter = -1;
      for (int f = 0; f < cell_num_faces_; ++f)
      {
        if (face_orientations[f] != FaceOrientation::OUTGOING) continue;

        // Set flags and counters
        out_face_counter++;
        const auto& face = cell_->faces_[f];
        const bool local = cell_transport_view_->IsFaceLocal(f);
        const bool boundary = not face.has_neighbor_;
        const int locality = cell_transport_view_->FaceLocality(f);

        if (not boundary and not local) ++deploc_face_counter;

        sweep_dependency_interface_.SetupOutgoingFace(
          f, cell_mapping_->NumFaceNodes(f), face.neighbor_id_, local, boundary, locality);

        aah_sweep_depinterf.out_face_counter = out_face_counter;
        aah_sweep_depinterf.deploc_face_counter = deploc_face_counter;

        OutgoingSurfaceOperations();
      } // for face

      ExecuteKernels(post_cell_dir_sweep_callbacks_);
    } // for n
  }   // for cell
}

const double*
AAH_SweepDependencyInterface::GetUpwindPsi(int face_node_local_idx) const
{
  const double* psi;
  if (on_local_face_)
    psi = fluds_->UpwindPsi(spls_index, in_face_counter, face_node_local_idx, 0, angle_set_index_);
  else if (not on_boundary_)
    psi = fluds_->NLUpwindPsi(preloc_face_counter, face_node_local_idx, 0, angle_set_index_);
  else
    psi = angle_set_->PsiBndry(neighbor_id_,
                               angle_num_,
                               cell_local_id_,
                               current_face_idx_,
                               face_node_local_idx,
                               gs_gi_,
                               gs_ss_begin_,
                               surface_source_active_);
  return psi;
}

double*
AAH_SweepDependencyInterface::GetDownwindPsi(int face_node_local_idx) const
{
  double* psi;
  if (on_local_face_)
    psi = fluds_->OutgoingPsi(spls_index, out_face_counter, face_node_local_idx, angle_set_index_);
  else if (not on_boundary_)
    psi = fluds_->NLOutgoingPsi(deploc_face_counter, face_node_local_idx, angle_set_index_);
  else if (is_reflecting_bndry_)
    psi = angle_set_->ReflectingPsiOutBoundBndry(neighbor_id_,
                                                 angle_num_,
                                                 cell_local_id_,
                                                 current_face_idx_,
                                                 face_node_local_idx,
                                                 gs_ss_begin_);
  else
    psi = nullptr;

  return psi;
}

} // namespace lbs
