#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_curvilinear_solver/sweep_chunks/lbs_curvilinear_sweep_chunk_pwl.h"

#include "framework/math/quadratures/curvilinear_angular_quadrature.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/logging/log_exceptions.h"

namespace opensn
{
namespace lbs
{

SweepChunkPWLRZ::SweepChunkPWLRZ(
  const MeshContinuum& grid,
  const SpatialDiscretization& discretization_primary,
  const std::vector<lbs::UnitCellMatrices>& unit_cell_matrices,
  const std::vector<lbs::UnitCellMatrices>& secondary_unit_cell_matrices,
  std::vector<lbs::CellLBSView>& cell_transport_views,
  std::vector<double>& destination_phi,
  std::vector<double>& destination_psi,
  const std::vector<double>& source_moments,
  lbs::LBSGroupset& groupset,
  const std::map<int, lbs::XSPtr>& xs,
  int num_moments,
  int max_num_cell_dofs)
  : AAH_SweepChunk(grid,
                   discretization_primary,
                   unit_cell_matrices,
                   cell_transport_views,
                   destination_phi,
                   destination_psi,
                   source_moments,
                   groupset,
                   xs,
                   num_moments,
                   max_num_cell_dofs),
    secondary_unit_cell_matrices_(secondary_unit_cell_matrices),
    unknown_manager_(),
    psi_sweep_(),
    normal_vector_boundary_()
{
  const auto curvilinear_product_quadrature =
    std::dynamic_pointer_cast<CurvilinearAngularQuadrature>(groupset_.quadrature_);

  if (!curvilinear_product_quadrature)
    throw std::invalid_argument("D_DO_RZ_SteadyState::SweepChunkPWL::SweepChunkPWL : "
                                "invalid angular quadrature");

  //  configure unknown manager for quantities that depend on polar level
  const size_t dir_map_size = curvilinear_product_quadrature->GetDirectionMap().size();
  for (size_t m = 0; m < dir_map_size; ++m)
    unknown_manager_.AddUnknown(UnknownType::VECTOR_N, groupset_.groups_.size());

  //  allocate storage for sweeping dependency
  const unsigned int n_dof = discretization_primary.GetNumLocalDOFs(unknown_manager_);
  psi_sweep_.resize(n_dof);

  //  initialise mappings from direction linear index
  for (const auto& dir_set : curvilinear_product_quadrature->GetDirectionMap())
    for (const auto& dir_idx : dir_set.second)
      map_polar_level_.emplace(dir_idx, dir_set.first);

  //  set normal vector for symmetric boundary condition
  const int d = (grid_.Attributes() & DIMENSION_1) ? 2 : 0;
  normal_vector_boundary_ = Vector3(0.0, 0.0, 0.0);
  normal_vector_boundary_(d) = 1;

  RegisterKernel("FEMRZVolumetricGradTerm",
                 std::bind(&SweepChunkPWLRZ::KernelFEMRZVolumetricGradientTerm, this));
  RegisterKernel("FEMRZUpwindSurfaceIntegrals",
                 std::bind(&SweepChunkPWLRZ::KernelFEMRZUpwindSurfaceIntegrals, this));

  // Setup callbacks
  cell_data_callbacks_.push_back(std::bind(&SweepChunkPWLRZ::CellDataCallback, this));

  direction_data_callbacks_and_kernels_ = {std::bind(&SweepChunkPWLRZ::DirectionDataCallback, this),
                                           Kernel("FEMRZVolumetricGradTerm")};

  surface_integral_kernels_ = {Kernel("FEMUpwindSurfaceIntegrals")};

  mass_term_kernels_ = {Kernel("FEMSSTDMassTerms")};

  // flux_update_kernels_ unchanged

  post_cell_dir_sweep_callbacks_.push_back(
    std::bind(&SweepChunkPWLRZ::PostCellDirSweepCallback, this));
}

void
SweepChunkPWLRZ::CellDataCallback()
{
  const auto& fe_intgrl_values_secondary = secondary_unit_cell_matrices_[cell_local_id_];

  Maux_ = &fe_intgrl_values_secondary.M_matrix;
}

void
SweepChunkPWLRZ::DirectionDataCallback()
{
  polar_level_ = map_polar_level_[direction_num_];
  const auto curvilinear_product_quadrature =
    std::dynamic_pointer_cast<CurvilinearAngularQuadrature>(groupset_.quadrature_);

  ChiLogicalErrorIf(not curvilinear_product_quadrature,
                    "Failure to cast angular quadrature to "
                    "CurvilinearAngularQuadrature");

  fac_diamond_difference_ =
    curvilinear_product_quadrature->GetDiamondDifferenceFactor()[direction_num_];
  fac_streaming_operator_ =
    curvilinear_product_quadrature->GetStreamingOperatorFactor()[direction_num_];
}

void
SweepChunkPWLRZ::PostCellDirSweepCallback()
{
  const auto f0 = 1 / fac_diamond_difference_;
  const auto f1 = f0 - 1;
  for (size_t i = 0; i < cell_num_nodes_; ++i)
  {
    const auto ir = grid_fe_view_.MapDOFLocal(*cell_, i, unknown_manager_, polar_level_, gs_gi_);
    for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
      psi_sweep_[ir + gsg] = f0 * b_[gsg][i] - f1 * psi_sweep_[ir + gsg];
  }
}

void
SweepChunkPWLRZ::KernelFEMRZVolumetricGradientTerm()
{
  const auto& G = *G_;
  const auto& Maux = *Maux_;

  for (int i = 0; i < cell_num_nodes_; ++i)
    for (int j = 0; j < cell_num_nodes_; ++j)
    {
      Amat_[i][j] = omega_.Dot(G[i][j]) + fac_streaming_operator_ * Maux[i][j];
      const auto jr = grid_fe_view_.MapDOFLocal(*cell_, j, unknown_manager_, polar_level_, gs_gi_);
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        b_[gsg][i] += fac_streaming_operator_ * Maux[i][j] * psi_sweep_[jr + gsg];
    }
}

void
SweepChunkPWLRZ::KernelFEMRZUpwindSurfaceIntegrals()
{
  const size_t f = sweep_dependency_interface_.current_face_idx_;
  if (sweep_dependency_interface_.on_boundary_)
  {
    const auto& face_normal = cell_->faces_[f].normal_;
    //------------------------------------------------------------------
    //  determine whether incoming direction is incident on the point
    //  of symmetry or on the axis of symmetry
    //  N.B.: a face is considered to be on the point/axis of symmetry
    //  if all are true:
    //    1. the face normal is antiparallel to $\vec{e}_{d}$
    //    2. all vertices of the face exhibit $v_{d} = 0$
    //  with $d = 2$ for 1D geometries and $d = 0$ for 2D geometries.
    //  Thanks to the verifications performed during initialisation,
    //  at this point it is necessary to confirm only the orientation.
    const bool incident_on_symmetric_boundary =
      (face_normal.Dot(normal_vector_boundary_) < -0.999999);
    if (incident_on_symmetric_boundary) return;
  }

  const auto& M_surf_f = (*M_surf_)[f];
  const double mu = face_mu_values_[f];
  const size_t num_face_nodes = cell_mapping_->NumFaceNodes(f);
  for (int fi = 0; fi < num_face_nodes; ++fi)
  {
    const int i = cell_mapping_->MapFaceNode(f, fi);
    for (int fj = 0; fj < num_face_nodes; ++fj)
    {
      const int j = cell_mapping_->MapFaceNode(f, fj);

      const double* psi = sweep_dependency_interface_.GetUpwindPsi(fj);

      const double mu_Nij = -mu * M_surf_f[i][j];
      Amat_[i][j] += mu_Nij;
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        b_[gsg][i] += psi[gsg] * mu_Nij;
    } // for face node j
  }   // for face node i
}

} // namespace lbs
} // namespace opensn
