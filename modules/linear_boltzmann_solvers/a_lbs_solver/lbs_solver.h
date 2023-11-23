#pragma once

#include "framework/physics/solver_base/solver.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_structs.h"
#include "framework/mesh/sweep_utilities/sweep_namespace.h"
#include "framework/mesh/sweep_utilities/sweep_boundary/sweep_boundary.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/point_source/point_source.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/distributed_source/distributed_source.h"

#include <petscksp.h>

namespace opensn
{

class MPICommunicatorSet;
class GridFaceHistogram;
typedef std::shared_ptr<GridFaceHistogram> GridFaceHistogramPtr;

class TimeIntegration;

namespace lbs
{
class AGSLinearSolver;
class WGSLinearSolver;
struct WGSContext;

/**Base class for all Linear Boltzmann Solvers.*/
class LBSSolver : public opensn::Solver
{
public:
  typedef std::shared_ptr<AGSLinearSolver> AGSLinSolverPtr;
  typedef std::shared_ptr<LinearSolver> LinSolvePtr;

public:
  explicit LBSSolver(const std::string& text_name);
  /**
   * Input parameters based construction.
   */
  explicit LBSSolver(const InputParameters& params);

  LBSSolver(const LBSSolver&) = delete;
  LBSSolver& operator=(const LBSSolver&) = delete;

  virtual ~LBSSolver() = default;

  /**
   * Returns the source event tag used for logging the time it takes to set source moments.
   */
  size_t GetSourceEventTag() const;

  /**
   * Returns the time at which the last restart was written.
   */
  double LastRestartWrite() const;

  /**
   * Returns a reference to the time at which the last restart was written.
   */
  double& LastRestartWrite();

  /**
   * Returns a reference to the solver options.
   */
  lbs::Options& Options();

  /**
   * Returns a constant reference to the solver options.
   */
  const lbs::Options& Options() const;

  static InputParameters OptionsBlock();
  static InputParameters BoundaryOptionsBlock();
  void SetOptions(const InputParameters& params);
  void SetBoundaryOptions(const InputParameters& params);

  /**
   * Returns the number of moments for the solver. This will only be non-zero after initialization.
   */
  size_t NumMoments() const;

  /**
   * Returns the number of groups for the solver. This will only be non-zero after initialization.
   */
  size_t NumGroups() const;

  /**
   * Returns the number of precursors for the solver. This will only be non-zero after
   * initialization.
   */
  size_t NumPrecursors() const;

  /**
   * Returns the maximum number of precursors defined on any material. This will only be
   * non-zero after initialization.
   */
  size_t MaxPrecursorsPerMaterial() const;

  /**
   * Adds a group to the list of groups. If group id < 0, the id will be logically derived from the
   * list size. If >= 0 the id will be set to the id specified.
   */
  void AddGroup(int id);
  const std::vector<LBSGroup>& Groups() const;

  /**
   * Adds a groupset to the list of groupsets. The groupset id will be logically
   * derived from the list size.
   */
  void AddGroupset();
  std::vector<LBSGroupset>& Groupsets();
  const std::vector<LBSGroupset>& Groupsets() const;

  /**
   * Adds a point source to the solver.
   */
  void AddPointSource(PointSource&& point_source);

  /**
   * Clears all the point sources from the solver.
   */
  void ClearPointSources();

  /**
   * Constant accessor to the list of point sources.
   */
  const std::vector<PointSource>& PointSources() const;

  /**
   * Adds a distributed source to the solver.
   */
  void AddDistributedSource(DistributedSource&& distributed_source);

  /**
   * Clears all the distributed sources from the solver.
   */
  void ClearDistributedSources();

  /**
   * Constant accessor to the list of distributed sources.
   */
  const std::vector<DistributedSource>& DistributedSources() const;

  /**
   * Returns a reference to the map of material ids to XSs.
   */
  const std::map<int, XSPtr>& GetMatID2XSMap() const;

  /**
   * Returns a reference to the map of material ids to Isotropic Srcs.
   */
  const std::map<int, IsotropicSrcPtr>& GetMatID2IsoSrcMap() const;

  /**
   * Obtains a reference to the grid.
   */
  const MeshContinuum& Grid() const;

  /**
   * Obtains a reference to the spatial discretization.
   */
  const class SpatialDiscretization& SpatialDiscretization() const;

  /**
   * Returns read-only access to the unit cell matrices.
   */
  const std::vector<UnitCellMatrices>& GetUnitCellMatrices() const;

  /**
   * Returns read-only access to the unit ghost cell matrices.
   */
  const std::map<uint64_t, UnitCellMatrices>& GetUnitGhostCellMatrices() const;

  /**
   * Returns a reference to the list of local cell transport views.
   */
  const std::vector<lbs::CellLBSView>& GetCellTransportViews() const;

  /**
   * Read/Write access to the boundary preferences.
   */
  std::map<uint64_t, BoundaryPreference>& BoundaryPreferences();

  /**
   * Obtains a reference to the unknown manager for flux-moments.
   */
  const class UnknownManager& UnknownManager() const;

  /**
   * Returns the local node count for the flux-moments data structures.
   */
  size_t LocalNodeCount() const;

  /**
   * Returns the global node count for the flux-moments data structures.
   */
  size_t GlobalNodeCount() const;

  /**
   * Read/write access to source moments vector.
   */
  std::vector<double>& QMomentsLocal();

  /**
   * Read access to source moments vector.
   */
  const std::vector<double>& QMomentsLocal() const;

  /**
   * Read/write access to exterior src moments vector.
   */
  std::vector<double>& ExtSrcMomentsLocal();

  /**
   * Read access to exterior src moments vector.
   */
  const std::vector<double>& ExtSrcMomentsLocal() const;

  /**
   * Read/write access to last updated flux vector.
   */
  std::vector<double>& PhiOldLocal();

  /**
   * Read access to last updated flux vector.
   */
  const std::vector<double>& PhiOldLocal() const;

  /**
   * Read/write access to newest updated flux vector.
   */
  std::vector<double>& PhiNewLocal();

  /**
   * Read access to newest updated flux vector.
   */
  const std::vector<double>& PhiNewLocal() const;

  /**
   * Read/write access to newest updated precursors vector.
   */
  std::vector<double>& PrecursorsNewLocal();

  /**
   * Read access to newest updated precursors vector.
   */
  const std::vector<double>& PrecursorsNewLocal() const;

  /**
   * Read/write access to newest updated angular flux vector.
   */
  std::vector<VecDbl>& PsiNewLocal();

  /**
   * Read access to newest updated angular flux vector.
   */
  const std::vector<VecDbl>& PsiNewLocal() const;

  /**
   * Returns the sweep boundaries as a read only reference
   */
  const std::map<uint64_t, std::shared_ptr<SweepBndry>>& SweepBoundaries() const;

  SetSourceFunction GetActiveSetSourceFunction() const;

  AGSLinSolverPtr GetPrimaryAGSSolver();

  std::vector<LinSolvePtr>& GetWGSSolvers();

  WGSContext& GetWGSContext(int groupset_id);

  /**
   * Gets the local and global number of iterative unknowns. This normally is only the flux moments,
   * however, the sweep based solvers might include delayed angular fluxes in this number.
   */
  virtual std::pair<size_t, size_t> GetNumPhiIterativeUnknowns();

  /**
   * Gets the local handle of a flux-moment based field function.
   */
  size_t MapPhiFieldFunction(size_t g, size_t m) const;

  /**
   * Returns the local handle to the power generation field function, if enabled.
   */
  size_t GetHandleToPowerGenFieldFunc() const;

  void Initialize() override;

  /**
   * Initializes default materials and physics materials.
   */
  void InitMaterials();

  /**Initializes the Within-Group DSA solver. */
  void InitWGDSA(LBSGroupset& groupset, bool vaccum_bcs_are_dirichlet = true);
  /**Creates a vector from a lbs primary stl vector where only the
   * scalar moments are mapped to the DOFs needed by WGDSA.*/
  std::vector<double> WGSCopyOnlyPhi0(const LBSGroupset& groupset,
                                      const std::vector<double>& phi_in);
  /**From the WGDSA DOFs, projects the scalar moments back into a
   * primary STL vector.*/
  void GSProjectBackPhi0(const LBSGroupset& groupset,
                         const std::vector<double>& input,
                         std::vector<double>& output);

  /**Assembles a delta-phi vector on the first moment.*/
  void AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in,
                                   std::vector<double>& delta_phi_local);

  /**
   * DAssembles a delta-phi vector on the first moment.
   */
  void DisAssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                      const std::vector<double>& delta_phi_local,
                                      std::vector<double>& ref_phi_new);

  /**
   * Assembles a delta-phi vector on the first moment.
   */
  void AssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                   const std::vector<double>& phi_in,
                                   std::vector<double>& delta_phi_local);
  /**
   * DAssembles a delta-phi vector on the first moment.
   */
  void DisAssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                      const std::vector<double>& delta_phi_local,
                                      std::vector<double>& ref_phi_new);

  /**
   * Writes phi_old to restart file.
   */
  void WriteRestartData(const std::string& folder_name, const std::string& file_base) const;

  /**
   * Read phi_old from restart file.
   */
  void ReadRestartData(const std::string& folder_name, const std::string& file_base);

  /**
   * Writes a full angular flux vector to file.
   */
  void WriteAngularFluxes(const std::vector<std::vector<double>>& src,
                          const std::string& file_base) const;

  /**
   * Reads a full angular flux vector from a file into the specified vector.
   */
  void ReadAngularFluxes(const std::string& file_base,
                         std::vector<std::vector<double>>& dest) const;

  /**
   * Writes a groupset angular flux vector to file.
   */
  void WriteGroupsetAngularFluxes(const LBSGroupset& groupset,
                                  const std::vector<double>& src,
                                  const std::string& file_base) const;

  /**
   * Reads the groupset angular fluxes from a file into the specified vector.
   */
  void ReadGroupsetAngularFluxes(const std::string& file_base,
                                 const LBSGroupset& groupset,
                                 std::vector<double>& dest) const;

  /**
   * Makes a source-moments vector from scattering and fission based on the latest phi-solution.
   */
  std::vector<double> MakeSourceMomentsFromPhi();

  /**
   * Writes a given flux moments vector to file.
   */
  void WriteFluxMoments(const std::vector<double>& src, const std::string& file_base) const;

  /**
   * Reads a flux moments vector from a file into the specified vector.
   */
  void ReadFluxMoments(const std::string& file_base,
                       std::vector<double>& dest,
                       bool single_file = false) const;

  /**
   * Copy relevant section of phi_old to the field functions.
   */
  void UpdateFieldFunctions();

  /**
   * Sets the internal phi vector to the value in the associated field function.
   */
  void SetPhiFromFieldFunctions(PhiSTLOption which_phi,
                                const std::vector<size_t>& m_indices,
                                const std::vector<size_t>& g_indices);

  /**
   * Compute the total fission production in the problem.
   * \author Zachary Hardy.
   */
  double ComputeFissionProduction(const std::vector<double>& phi);

  /**
   * Computes the total fission rate in the problem.
   * \author Zachary Hardy.
   */
  double ComputeFissionRate(const std::vector<double>& phi);

  /**
   * Compute the steady state delayed neutron precursor concentrations.
   */
  void ComputePrecursors();

  /**
   * Sets a value to the zeroth (scalar) moment of the vector.
   */
  virtual void SetPhiVectorScalarValues(std::vector<double>& phi_vector, double value);
  /**
   * Scales a flux moment vector. For sweep methods the delayed angular fluxes will also be scaled.
   */
  virtual void ScalePhiVector(PhiSTLOption which_phi, double value);
  /**
   * Assembles a vector for a given groupset from a source vector.
   */
  virtual void
  SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x, PhiSTLOption which_phi);

  /**
   * Assembles a vector for a given groupset from a source vector.
   */
  virtual void
  SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset, Vec x_src, PhiSTLOption which_phi);

  /**
   * Assembles a vector for a given groupset from a source vector.
   */
  virtual void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                             const std::vector<double>& x_src,
                                             std::vector<double>& y);

  /**
   * Assembles a vector for a given groupset from a source vector.
   */
  virtual void GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                             PhiSTLOption from_which_phi,
                                             PhiSTLOption to_which_phi);

  /**
   * Assembles a vector for a given group span from a source vector.
   */
  virtual void SetGroupScopedPETScVecFromPrimarySTLvector(int first_group_id,
                                                          int last_group_id,
                                                          Vec x,
                                                          const std::vector<double>& y);

  /**
   * Assembles a vector for a given groupset from a source vector.
   */
  virtual void SetPrimarySTLvectorFromGroupScopedPETScVec(int first_group_id,
                                                          int last_group_id,
                                                          Vec x_src,
                                                          std::vector<double>& y);

  /**
   * Assembles a PETSc vector from multiple groupsets.
   */
  virtual void SetMultiGSPETScVecFromPrimarySTLvector(const std::vector<int>& gs_ids,
                                                      Vec x,
                                                      PhiSTLOption which_phi);

  /**
   * Disassembles a multiple Groupset PETSc vector STL vectors.
   */
  virtual void SetPrimarySTLvectorFromMultiGSPETScVecFrom(const std::vector<int>& gs_ids,
                                                          Vec x_src,
                                                          PhiSTLOption which_phi);

protected:
  /**
   * Performs general input checks before initialization continues.
   */
  virtual void PerformInputChecks();

  /**
   * Prints header information of simulation.
   */
  void PrintSimHeader();

  virtual void InitializeSpatialDiscretization();
  void ComputeUnitIntegrals();

  /**
   * Initializes common groupset items.
   */
  void InitializeGroupsets();

  /**
   * Computes the number of moments for the given mesher types
   */
  void ComputeNumberOfMoments();

  /**
   * Initializes parallel arrays.
   */
  virtual void InitializeParrays();

  void InitializeFieldFunctions();

  /**Initializes transport related boundaries. */
  void InitializeBoundaries();
  virtual void InitializeSolverSchemes();
  virtual void InitializeWGSSolvers(){};
  /**Initializes the Within-Group DSA solver. */
  void InitTGDSA(LBSGroupset& groupset);

  size_t source_event_tag_ = 0;
  double last_restart_write_ = 0.0;

  lbs::Options options_;
  size_t num_moments_ = 0;
  size_t num_groups_ = 0;
  size_t num_precursors_ = 0;
  size_t max_precursors_per_material_ = 0;

  std::vector<LBSGroup> groups_;
  std::vector<LBSGroupset> groupsets_;

  std::map<int, XSPtr> matid_to_xs_map_;
  std::map<int, IsotropicSrcPtr> matid_to_src_map_;

  std::vector<PointSource> point_sources_;
  std::vector<DistributedSource> distributed_sources_;

  std::shared_ptr<MeshContinuum> grid_ptr_;
  std::shared_ptr<opensn::SpatialDiscretization> discretization_ = nullptr;

  std::vector<CellFaceNodalMapping> grid_nodal_mappings_;
  std::shared_ptr<MPICommunicatorSet> grid_local_comm_set_ = nullptr;
  GridFaceHistogramPtr grid_face_histogram_ = nullptr;

  std::vector<UnitCellMatrices> unit_cell_matrices_;
  std::map<uint64_t, UnitCellMatrices> unit_ghost_cell_matrices_;
  std::vector<lbs::CellLBSView> cell_transport_views_;

  std::map<uint64_t, BoundaryPreference> boundary_preferences_;
  std::map<uint64_t, std::shared_ptr<SweepBndry>> sweep_boundaries_;

  opensn::UnknownManager flux_moments_uk_man_;

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
  std::shared_ptr<const TimeIntegration> time_integration_ = nullptr;

  /**Cleans up memory consuming items. */
  static void CleanUpWGDSA(LBSGroupset& groupset);
  /**Cleans up memory consuming items. */
  static void CleanUpTGDSA(LBSGroupset& groupset);

public:
  static std::map<std::string, uint64_t> supported_boundary_names;
  static std::map<uint64_t, std::string> supported_boundary_ids;

  /**
   * Returns the input parameters for this object.
   */
  static InputParameters GetInputParameters();
};

} // namespace lbs
} // namespace opensn
