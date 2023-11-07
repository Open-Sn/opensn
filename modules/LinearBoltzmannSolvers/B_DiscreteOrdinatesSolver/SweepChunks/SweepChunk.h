#pragma once

#include "framework/mesh/sweep_utilities/sweep_chunk_base.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_structs.h"

namespace chi_math
{
class CellMapping;
}

namespace lbs
{

struct SweepDependencyInterface
{
  size_t groupset_angle_group_stride_;
  size_t groupset_group_stride_;

  chi_mesh::sweep_management::AngleSet* angle_set_ = nullptr;
  bool surface_source_active_ = false;

  size_t gs_ss_begin_ = 0;
  int gs_gi_ = 0;

  const chi_mesh::Cell* cell_ptr_ = nullptr;
  uint64_t cell_local_id_ = 0;

  size_t angle_set_index_ = 0;
  size_t angle_num_ = 0;

public: // Set using SetupIncomingFace
  int current_face_idx_ = 0;
  size_t num_face_nodes_ = 0;
  uint64_t neighbor_id_ = 0;
  int face_locality_ = 0;

  bool on_local_face_ = false;
  bool on_boundary_ = false;

public:
  bool is_reflecting_bndry_ = false;

  SweepDependencyInterface() = default;

  virtual const double* GetUpwindPsi(int face_node_local_idx) const = 0;
  virtual double* GetDownwindPsi(int face_node_local_idx) const = 0;

  /**Sets data for the current incoming face.*/
  virtual void SetupIncomingFace(
    int face_id, size_t num_face_nodes, uint64_t neighbor_id, bool on_local_face, bool on_boundary);
  /**Sets data for the current outgoing face.*/
  virtual void SetupOutgoingFace(int face_id,
                                 size_t num_face_nodes,
                                 uint64_t neighbor_id,
                                 bool on_local_face,
                                 bool on_boundary,
                                 int locality);

  virtual ~SweepDependencyInterface() = default;
};

/**Base class for LBS sweepers*/
class SweepChunk : public chi_mesh::sweep_management::SweepChunk
{
public:
  SweepChunk(std::vector<double>& destination_phi,
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
             std::unique_ptr<SweepDependencyInterface> sweep_dependency_interface_ptr);

protected:
  typedef std::function<void()> CallbackFunction;

  const chi_mesh::MeshContinuum& grid_;
  const chi_math::SpatialDiscretization& grid_fe_view_;
  const std::vector<UnitCellMatrices>& unit_cell_matrices_;
  std::vector<lbs::CellLBSView>& grid_transport_view_;
  const std::vector<double>& q_moments_;
  const LBSGroupset& groupset_;
  const std::map<int, XSPtr>& xs_;
  const int num_moments_;
  const bool save_angular_flux_;

  std::unique_ptr<SweepDependencyInterface> sweep_dependency_interface_ptr_;
  SweepDependencyInterface& sweep_dependency_interface_;

  const size_t groupset_angle_group_stride_;
  const size_t groupset_group_stride_;

  // Runtime params
  size_t gs_ss_size_ = 0;
  size_t gs_ss_begin_ = 0;
  int gs_gi_ = 0;

  std::vector<std::vector<double>> Amat_;
  std::vector<std::vector<double>> Atemp_;
  std::vector<double> source_;
  std::vector<std::vector<double>> b_;

  // Cell items
  uint64_t cell_local_id_ = 0;
  const chi_mesh::Cell* cell_ = nullptr;
  const chi_math::CellMapping* cell_mapping_ = nullptr;
  CellLBSView* cell_transport_view_ = nullptr;
  size_t cell_num_faces_ = 0;
  size_t cell_num_nodes_ = 0;
  const MatVec3* G_ = nullptr;
  const MatDbl* M_ = nullptr;
  const std::vector<MatDbl>* M_surf_ = nullptr;
  const std::vector<VecDbl>* IntS_shapeI_ = nullptr;

  /**Callbacks at phase 1 : cell data established*/
  std::vector<CallbackFunction> cell_data_callbacks_;

  std::vector<double> face_mu_values_;
  size_t direction_num_ = 0;
  chi_mesh::Vector3 omega_;
  double direction_qweight_ = 0.0;

  /**Callbacks at phase 2 : direction data established*/
  std::vector<CallbackFunction> direction_data_callbacks_and_kernels_;

  /**Callbacks at phase 3 : Surface integrals*/
  std::vector<CallbackFunction> surface_integral_kernels_;

  size_t g_ = 0;
  size_t gsg_ = 0;
  double sigma_tg_ = 0.0;

  /**Callbacks at phase 4 : group by group mass terms*/
  std::vector<CallbackFunction> mass_term_kernels_;

  /**Callbacks at phase 5 : flux updates*/
  std::vector<CallbackFunction> flux_update_kernels_;

  /**Callbacks at phase 6 : Post cell-dir sweep*/
  std::vector<CallbackFunction> post_cell_dir_sweep_callbacks_;

  // 02 operations
  /**Registers a kernel as a named callback function*/
  void RegisterKernel(const std::string& name, CallbackFunction function);
  /**Returns a kernel if the given name exists.*/
  CallbackFunction Kernel(const std::string& name) const;
  /**Executes the supplied kernels list.*/
  static void ExecuteKernels(const std::vector<CallbackFunction>& kernels);
  /**Operations when outgoing fluxes are handled including passing
   * face angular fluxes downstream and computing
   * balance parameters (i.e. outflow)
   * */
  virtual void OutgoingSurfaceOperations();

  // kernels
public:
  /**Assembles the volumetric gradient term.*/
  void KernelFEMVolumetricGradientTerm();
  /**Performs the integral over the surface of a face.*/
  void KernelFEMUpwindSurfaceIntegrals();
  /**Assembles angular sources and applies the mass matrix terms.*/
  void KernelFEMSTDMassTerms();
  /**Adds a single direction's contribution to the moment integrals.*/
  void KernelPhiUpdate();
  /**Updates angular fluxes.*/
  void KernelPsiUpdate();

private:
  std::map<std::string, CallbackFunction> kernels_;
};

} // namespace lbs
