#pragma once

#include "physics/SolverBase/chi_solver.h"

#include "A_LBSSolver/Groupset/lbs_groupset.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/LinearSolver/linear_solver.h"
#include "lbs_structs.h"
#include "mesh/SweepUtilities/sweep_namespace.h"
#include "mesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"

#include "A_LBSSolver/PointSource/lbs_point_source.h"

#include <petscksp.h>

namespace lbs
{
template <class MatType, class VecType, class SolverType>
class AGSLinearSolver;
template <class MatType, class VecType, class SolverType>
class WGSLinearSolver;
template <class MatType, class VecType, class SolverType>
struct WGSContext;
} // namespace lbs

namespace chi
{
class ChiMPICommunicatorSet;
}
typedef std::shared_ptr<chi::ChiMPICommunicatorSet> MPILocalCommSetPtr;

namespace chi_mesh
{
class GridFaceHistogram;
}
typedef std::shared_ptr<chi_mesh::GridFaceHistogram> GridFaceHistogramPtr;

namespace chi_math
{
class TimeIntegration;
}

namespace lbs
{

// ################################################################### Class def
/**Base class for all Linear Boltzmann Solvers.*/
class LBSSolver : public chi_physics::Solver
{
public:
  typedef std::shared_ptr<AGSLinearSolver<Mat, Vec, KSP>> AGSLinSolverPtr;
  typedef std::shared_ptr<chi_math::LinearSolver<Mat, Vec, KSP>> LinSolvePtr;

protected:
  typedef chi_mesh::sweep_management::CellFaceNodalMapping CellFaceNodalMapping;

  size_t source_event_tag_ = 0;
  double last_restart_write_ = 0.0;

  lbs::Options options_;
  size_t num_moments_ = 0;
  size_t num_groups_ = 0;
  size_t num_precursors_ = 0;
  size_t max_precursors_per_material_ = 0;

  std::vector<LBSGroup> groups_;
  std::vector<LBSGroupset> groupsets_;
  std::vector<PointSource> point_sources_;

  std::map<int, XSPtr> matid_to_xs_map_;
  std::map<int, IsotropicSrcPtr> matid_to_src_map_;

  std::shared_ptr<chi_math::SpatialDiscretization> discretization_ = nullptr;
  chi_mesh::MeshContinuumPtr grid_ptr_;

  std::vector<CellFaceNodalMapping> grid_nodal_mappings_;
  MPILocalCommSetPtr grid_local_comm_set_ = nullptr;
  GridFaceHistogramPtr grid_face_histogram_ = nullptr;

  std::vector<UnitCellMatrices> unit_cell_matrices_;
  std::map<uint64_t, UnitCellMatrices> unit_ghost_cell_matrices_;
  std::vector<lbs::CellLBSView> cell_transport_views_;

  std::map<uint64_t, BoundaryPreference> boundary_preferences_;
  std::map<uint64_t, std::shared_ptr<SweepBndry>> sweep_boundaries_;

  chi_math::UnknownManager flux_moments_uk_man_;

  size_t max_cell_dof_count_ = 0;
  uint64_t local_node_count_ = 0;
  uint64_t glob_node_count_ = 0;

  std::vector<double> q_moments_local_, ext_src_moments_local_;
  std::vector<double> phi_new_local_, phi_old_local_;
  std::vector<std::vector<double>> psi_new_local_;
  std::vector<double> precursor_new_local_;

  SetSourceFunction active_set_source_function_;

  std::vector<AGSLinSolverPtr> ags_solvers_;
  std::vector<LinSolvePtr> wgs_solvers_;
  AGSLinSolverPtr primary_ags_solver_;

  std::map<std::pair<size_t, size_t>, size_t> phi_field_functions_local_map_;
  size_t power_gen_fieldfunc_local_handle_ = 0;

  /**Time integration parameter meant to be set by an executor*/
  std::shared_ptr<const chi_math::TimeIntegration> time_integration_ = nullptr;

public:
  static chi::InputParameters GetInputParameters();
  explicit LBSSolver(const std::string& text_name);
  explicit LBSSolver(const chi::InputParameters& params);

  LBSSolver(const LBSSolver&) = delete;
  LBSSolver& operator=(const LBSSolver&) = delete;

  virtual ~LBSSolver() = default;

  size_t GetSourceEventTag() const;

  double LastRestartWrite() const;
  double& LastRestartWrite();

  lbs::Options& Options();
  const lbs::Options& Options() const;

  static chi::InputParameters OptionsBlock();
  static chi::InputParameters BoundaryOptionsBlock();
  void SetOptions(const chi::InputParameters& params);
  void SetBoundaryOptions(const chi::InputParameters& params);

  size_t NumMoments() const;
  size_t NumGroups() const;
  size_t NumPrecursors() const;

  size_t GetMaxPrecursorsPerMaterial() const;

  void AddGroup(int id);
  const std::vector<LBSGroup>& Groups() const;

  void AddGroupset();
  std::vector<LBSGroupset>& Groupsets();
  const std::vector<LBSGroupset>& Groupsets() const;

  void AddPointSource(PointSource psrc);
  void ClearPointSources();
  const std::vector<PointSource>& PointSources() const;

  const std::map<int, XSPtr>& GetMatID2XSMap() const;
  const std::map<int, IsotropicSrcPtr>& GetMatID2IsoSrcMap() const;

  const chi_math::SpatialDiscretization& SpatialDiscretization() const;
  const std::vector<UnitCellMatrices>& GetUnitCellMatrices() const;
  const chi_mesh::MeshContinuum& Grid() const;

  const std::vector<lbs::CellLBSView>& GetCellTransportViews() const;

  std::map<uint64_t, BoundaryPreference>& BoundaryPreferences();

  const chi_math::UnknownManager& UnknownManager() const;

  size_t LocalNodeCount() const;
  size_t GlobalNodeCount() const;

  std::vector<double>& QMomentsLocal();
  const std::vector<double>& QMomentsLocal() const;
  std::vector<double>& ExtSrcMomentsLocal();
  const std::vector<double>& ExtSrcMomentsLocal() const;
  std::vector<double>& PhiOldLocal();
  const std::vector<double>& PhiOldLocal() const;
  std::vector<double>& PhiNewLocal();
  const std::vector<double>& PhiNewLocal() const;
  std::vector<double>& PrecursorsNewLocal();
  const std::vector<double>& PrecursorsNewLocal() const;
  std::vector<VecDbl>& PsiNewLocal();
  const std::vector<VecDbl>& PsiNewLocal() const;

  /**Returns the sweep boundaries as a read only reference*/
  const std::map<uint64_t, std::shared_ptr<SweepBndry>>& SweepBoundaries() const;

  SetSourceFunction GetActiveSetSourceFunction() const;

  AGSLinSolverPtr GetPrimaryAGSSolver();

  std::vector<LinSolvePtr>& GetWGSSolvers();

  WGSContext<Mat, Vec, KSP>& GetWGSContext(int groupset_id);

  virtual std::pair<size_t, size_t> GetNumPhiIterativeUnknowns();

  /**Gets the local handle of a flux-moment based field function.*/
  size_t MapPhiFieldFunction(size_t g, size_t m) const;

  /**Returns the local handle to the power generation field function, if
   * enabled.*/
  size_t GetHandleToPowerGenFieldFunc() const;

  // 01
  void Initialize() override;

protected:
  // 01a
  virtual void PerformInputChecks();
  // 01b
  void PrintSimHeader();

public:
  // 01c
  void InitMaterials();

protected:
  // 01d
  virtual void InitializeSpatialDiscretization();
  void ComputeUnitIntegrals();
  // 01e
  void InitializeGroupsets();
  // 01f
  void ComputeNumberOfMoments();
  // 01g
  virtual void InitializeParrays();
  //   a
  void InitializeFieldFunctions();
  // 01h
  void InitializeBoundaries();

public:
  // 01i
  void InitializePointSources();

protected:
  // 01j
  virtual void InitializeSolverSchemes();
  virtual void InitializeWGSSolvers(){};

  // 03d
public:
  void InitWGDSA(LBSGroupset& groupset, bool vaccum_bcs_are_dirichlet = true);
  std::vector<double> WGSCopyOnlyPhi0(const LBSGroupset& groupset,
                                      const std::vector<double>& phi_in);
  void GSProjectBackPhi0(const LBSGroupset& groupset,
                         const std::vector<double>& input,
                         std::vector<double>& output);
  //  std::vector<double> WGDSAMake
public:
  void AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in,
                                   std::vector<double>& delta_phi_local);

  void DisAssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                      const std::vector<double>& delta_phi_local,
                                      std::vector<double>& ref_phi_new);

protected:
  static void CleanUpWGDSA(LBSGroupset& groupset);

  // 03e
  void InitTGDSA(LBSGroupset& groupset);

public:
  void AssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in,
                                   std::vector<double>& delta_phi_local);
  void DisAssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                      const std::vector<double>& delta_phi_local,
                                      std::vector<double>& ref_phi_new);

protected:
  static void CleanUpTGDSA(LBSGroupset& groupset);
  // 04 File IO
public:
  // 04a
  void WriteRestartData(const std::string& folder_name, const std::string& file_base);
  void ReadRestartData(const std::string& folder_name, const std::string& file_base);
  // 04b
  void WriteGroupsetAngularFluxes(const LBSGroupset& groupset, const std::string& file_base);
  void ReadGroupsetAngularFluxes(LBSGroupset& groupset, const std::string& file_base);

  // 04c
  std::vector<double> MakeSourceMomentsFromPhi();
  void WriteFluxMoments(const std::string& file_base, const std::vector<double>& flux_moments);
  void ReadFluxMoments(const std::string& file_base,
                       std::vector<double>& flux_moments,
                       bool single_file = false);

  // 05a
  void UpdateFieldFunctions();
  void SetPhiFromFieldFunctions(PhiSTLOption which_phi,
                                const std::vector<size_t>& m_indices,
                                const std::vector<size_t>& g_indices);

  // 06b
public:
  double ComputeFissionProduction(const std::vector<double>& phi);
  double ComputeFissionRate(const std::vector<double>& phi);

  // 06c
public:
  void ComputePrecursors();

  // 07 Vector assembly
public:
  virtual void SetPhiVectorScalarValues(std::vector<double>& phi_vector, double value);
  virtual void ScalePhiVector(PhiSTLOption which_phi, double value);
  virtual void
  SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x, PhiSTLOption which_phi);

  virtual void
  SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset, Vec x_src, PhiSTLOption which_phi);

  virtual void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                             const std::vector<double>& x_src,
                                             std::vector<double>& y);

  virtual void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                             PhiSTLOption from_which_phi,
                                             PhiSTLOption to_which_phi);

  virtual void SetGroupScopedPETScVecFromPrimarySTLvector(int first_group_id,
                                                          int last_group_id,
                                                          Vec x,
                                                          const std::vector<double>& y);

  virtual void SetPrimarySTLvectorFromGroupScopedPETScVec(int first_group_id,
                                                          int last_group_id,
                                                          Vec x_src,
                                                          std::vector<double>& y);

  virtual void SetMultiGSPETScVecFromPrimarySTLvector(const std::vector<int>& gs_ids,
                                                      Vec x,
                                                      PhiSTLOption which_phi);

  virtual void SetPrimarySTLvectorFromMultiGSPETScVecFrom(const std::vector<int>& gs_ids,
                                                          Vec x_src,
                                                          PhiSTLOption which_phi);
};

} // namespace lbs


