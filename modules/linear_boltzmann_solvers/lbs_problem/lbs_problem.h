// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/source_functions/source_flags.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/point_source/point_source.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/volumetric_source/volumetric_source.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_view.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/math/geometry.h"
#include "framework/utils/hdf_utils.h"
#include <memory>
#include <petscksp.h>
#include <chrono>
#include <functional>

namespace opensn
{

class MPICommunicatorSet;
class GridFaceHistogram;
class FieldFunctionGridBased;
class TotalXSCarrier;
class OutflowCarrier;
class MeshCarrier;
class LBSSolverIO;
template <typename T>
class MemoryPinner;

/// Base class for all Linear Boltzmann Solvers.
class LBSProblem : public Problem, public std::enable_shared_from_this<LBSProblem>
{
public:
  using RestartDataHook = std::function<bool(hid_t)>;

  friend class LBSSolverIO;

  LBSProblem(const LBSProblem&) = delete;

  LBSProblem& operator=(const LBSProblem&) = delete;

  ~LBSProblem() override;

  /// Returns a constant reference to the solver options.
  const LBSOptions& GetOptions() const;

  /// Returns simulation time in seconds for time dependent problems.
  double GetTime() const;

  /// Sets simulation time in seconds for time dependent problems.
  void SetTime(double time);

  /// Sets dt.
  void SetTimeStep(double dt);

  /// Returns dt.
  double GetTimeStep() const;

  /// Sets theta for time discretization.
  void SetTheta(double theta);

  /// Returns theta for time discretization.
  double GetTheta() const;

  /// Returns true if the problem is currently in time-dependent mode.
  virtual bool IsTimeDependent() const;

  /// Set the problem to time-dependent mode.
  virtual void SetTimeDependentMode();

  /// Set the problem to steady-state mode.
  virtual void SetSteadyStateMode();

  /**
   * Toggle forward/adjoint transport mode.
   *
   * If the requested mode differs from the current mode, this performs a
   * mode-transition reset. Materials are reinitialized in the selected mode,
   * sources and boundaries are cleared, and solution vectors are zeroed.
   */
  void SetAdjoint(bool adjoint = true);

  /// Set the problem to forward mode.
  void SetForward();

  /// Returns true if the problem is in adjoint mode.
  bool IsAdjoint() const;

  void ZeroPhi();
  void CopyPhiNewToOld();
  void SetPhiOldFrom(const std::vector<double>& phi_old);
  void SetPhiNewFrom(const std::vector<double>& phi_new);
  void ScalePhiOld(double factor);
  void ScalePhiNew(double factor);
  void ZeroQMoments();
  void ScaleQMoments(double factor);
  void SetQMomentsFrom(const std::vector<double>& q_moments);
  void ScalePrecursors(double factor);
  void ZeroPrecursors();
  void ZeroExtSrcMoments();
  void ScaleExtSrcMoments(double factor);
  GeometryType GetGeometryType() const;

  /// Returns the number of moments for the solver.
  unsigned int GetNumMoments() const;

  unsigned int GetMaxCellDOFCount() const;

  unsigned int GetMinCellDOFCount() const;

  bool UseGPUs() const;

  /// Returns the number of groups for the solver.
  unsigned int GetNumGroups() const;

  /// Returns the scattering order for the solver.
  unsigned int GetScatteringOrder() const;

  /// Returns the number of precursors for the solver.
  unsigned int GetNumPrecursors() const;

  /// Returns the maximum number of precursors defined on any material.
  unsigned int GetMaxPrecursorsPerMaterial() const;

  const std::vector<LBSGroupset>& GetGroupsets() const;
  LBSGroupset& GetGroupset(size_t groupset_id);
  const LBSGroupset& GetGroupset(size_t groupset_id) const;
  size_t GetNumGroupsets() const;

  /// Adds a point source to the solver.
  void AddPointSource(std::shared_ptr<PointSource> point_source);

  /// Clears all the point sources from the solver.
  void ClearPointSources();

  /// Constant accessor to the list of point sources.
  const std::vector<std::shared_ptr<PointSource>>& GetPointSources() const;

  /// Adds a volumetric source to the solver.
  void AddVolumetricSource(std::shared_ptr<VolumetricSource> volumetric_source);

  /// Clears all the volumetric sources from the solver.
  void ClearVolumetricSources();

  /// Constant accessor to the list of volumetric sources.
  const std::vector<std::shared_ptr<VolumetricSource>>& GetVolumetricSources() const;

  /// Clears all the boundary conditions from the solver.
  virtual void ClearBoundaries() = 0;

  /// Returns a reference to the map of material ids to XSs.
  const BlockID2XSMap& GetBlockID2XSMap() const;

  /// Replaces the map of block ids to XSs and refreshes material data.
  virtual void SetBlockID2XSMap(const BlockID2XSMap& xs_map);

  /// Obtains a reference to the grid.
  std::shared_ptr<MeshContinuum> GetGrid() const;

  TotalXSCarrier* GetTotalXSCarrier() { return total_xs_carrier_.get(); }
  const TotalXSCarrier* GetTotalXSCarrier() const { return total_xs_carrier_.get(); }

  OutflowCarrier* GetOutflowCarrier() { return outflow_carrier_.get(); }
  const OutflowCarrier* GetOutflowCarrier() const { return outflow_carrier_.get(); }

  MeshCarrier* GetMeshCarrier() { return mesh_carrier_.get(); }
  const MeshCarrier* GetMeshCarrier() const { return mesh_carrier_.get(); }

  MemoryPinner<double>* GetSourceMomentsPinner() { return source_pinner_.get(); }
  const MemoryPinner<double>* GetSourceMomentsPinner() const { return source_pinner_.get(); }

  MemoryPinner<double>* GetPhiPinner() { return phi_pinner_.get(); }
  const MemoryPinner<double>* GetPhiPinner() const { return phi_pinner_.get(); }

  /// Obtains a reference to the spatial discretization.
  const SpatialDiscretization& GetSpatialDiscretization() const;

  /// Returns read-only access to the unit cell matrices.
  const std::vector<UnitCellMatrices>& GetUnitCellMatrices() const;

  /// Returns read-only access to the unit ghost cell matrices.
  const std::map<uint64_t, UnitCellMatrices>& GetUnitGhostCellMatrices() const;

  /// Returns a reference to the list of local cell transport views.
  std::vector<CellLBSView>& GetCellTransportViews();

  /// Returns a const reference to the list of local cell transport views.
  const std::vector<CellLBSView>& GetCellTransportViews() const;

  /// Obtains a reference to the unknown manager for flux-moments.
  const UnknownManager& GetUnknownManager() const;

  /// Returns the local node count for the flux-moments data structures.
  size_t GetLocalNodeCount() const;

  /// Returns the global node count for the flux-moments data structures.
  std::uint64_t GetGlobalNodeCount() const;

  /// Read/write access to source moments vector.
  std::vector<double>& GetQMomentsLocal();

  /// Read access to source moments vector.
  const std::vector<double>& GetQMomentsLocal() const;

  /// Read/write access to exterior src moments vector.
  std::vector<double>& GetExtSrcMomentsLocal();

  /// Read access to exterior src moments vector.
  const std::vector<double>& GetExtSrcMomentsLocal() const;

  /// Replaces exterior source moments vector.
  void SetExtSrcMomentsFrom(const std::vector<double>& ext_src_moments);

  /// Read/write access to last updated flux vector.
  std::vector<double>& GetPhiOldLocal();

  /// Read access to last updated flux vector.
  const std::vector<double>& GetPhiOldLocal() const;

  /// Read/write access to newest updated flux vector.
  std::vector<double>& GetPhiNewLocal();

  /// Read access to newest updated flux vector.
  const std::vector<double>& GetPhiNewLocal() const;

  /// Read/write access to newest updated precursors vector.
  std::vector<double>& GetPrecursorsNewLocal();

  /// Read access to newest updated precursors vector.
  const std::vector<double>& GetPrecursorsNewLocal() const;

  SetSourceFunction GetActiveSetSourceFunction() const;

  /**
   * Gets the local and global number of iterative unknowns. This normally is only the flux moments,
   * however, the sweep based solvers might include delayed angular fluxes in this number.
   */
  virtual std::pair<std::uint64_t, std::uint64_t> GetNumPhiIterativeUnknowns();

  /**
   * Returns a flux-moment field function for a given group and moment.
   *
   * The returned field function is created on demand from the current \c phi_new state.
   */
  std::shared_ptr<FieldFunctionGridBased> CreateScalarFluxFieldFunction(unsigned int g,
                                                                        unsigned int m = 0);

  /// Creates a named field function from a 1D XS or the special case `power`.
  std::shared_ptr<FieldFunctionGridBased> CreateFieldFunction(
    const std::string& name, const std::string& xs_name, double power_normalization_target = -1.0);

  bool ReadRestartData(const RestartDataHook& extra_reader = {});
  bool WriteRestartData(const RestartDataHook& extra_writer = {});

  bool TriggerRestartDump() const
  {
    if (options_.write_restart_time_interval <= std::chrono::seconds(0))
      return false;

    auto elapsed = std::chrono::system_clock::now() - options_.last_restart_write_time;
    return elapsed >= options_.write_restart_time_interval;
  }

  void UpdateRestartWriteTime()
  {
    options_.last_restart_write_time = std::chrono::system_clock::now();
  }

  /// Makes a source-moments vector from scattering and fission based on the latest phi-solution.
  std::vector<double> MakeSourceMomentsFromPhi();

  /**
   * A method for post-processing an adjoint solution.
   *
   * @note This does nothing for diffusion-based solvers.
   */
  virtual void ReorientAdjointSolution() {};

protected:
  /// Input parameters based construction.
  explicit LBSProblem(const InputParameters& params);

  /// Build runtime data structures once all constructor-time configuration is complete.
  void BuildRuntime();

  virtual void PrintSimHeader();

  virtual void ResetDerivedSolutionVectors() {}

  void ComputeUnitIntegrals();

  virtual void InitializeSpatialDiscretization();

  /// Initializes boundaries.
  virtual void InitializeBoundaries() {}

  /// Derived problems handle boundary options.
  virtual void SetBoundaryOptions(const InputParameters& params) = 0;

  void SetActiveSetSourceFunction(SetSourceFunction source_function);

  std::shared_ptr<FieldFunctionGridBased> CreateEmptyFieldFunction(const std::string& name) const;
  std::string MakeFieldFunctionName(const std::string& base_name) const;
  void UpdateScalarFluxFieldFunction(FieldFunctionGridBased& ff, unsigned int g, unsigned int m);
  void UpdateDerivedFieldFunction(FieldFunctionGridBased& ff,
                                  const std::string& xs_name,
                                  double power_normalization_target);
  virtual bool ReadProblemRestartData(hid_t file_id);
  virtual bool WriteProblemRestartData(hid_t file_id) const;

  LBSOptions options_;
  double time_ = 0.0;
  double theta_ = 1.0;
  double dt_ = 1.0;
  GeometryType geometry_type_ = GeometryType::INVALID;
  unsigned int num_moments_ = 0;
  unsigned int num_groups_ = 0;
  unsigned int scattering_order_ = 0;
  unsigned int num_precursors_ = 0;
  unsigned int max_precursors_per_material_ = 0;

  std::vector<LBSGroupset> groupsets_;

  BlockID2XSMap block_id_to_xs_map_;

  std::vector<std::shared_ptr<PointSource>> point_sources_;
  std::vector<std::shared_ptr<VolumetricSource>> volumetric_sources_;

  std::shared_ptr<MeshContinuum> grid_;
  std::shared_ptr<SpatialDiscretization> discretization_ = nullptr;

  std::vector<CellFaceNodalMapping> grid_nodal_mappings_;
  std::shared_ptr<MPICommunicatorSet> grid_local_comm_set_ = nullptr;

  std::vector<UnitCellMatrices> unit_cell_matrices_;
  std::map<uint64_t, UnitCellMatrices> unit_ghost_cell_matrices_;
  std::vector<CellLBSView> cell_transport_views_;

  UnknownManager flux_moments_uk_man_;

  unsigned int max_cell_dof_count_ = 0;
  unsigned int min_cell_dof_count_ = 0;
  uint64_t local_node_count_ = 0;
  std::uint64_t global_node_count_ = 0;

  std::vector<double> q_moments_local_, ext_src_moments_local_;
  std::vector<double> phi_new_local_, phi_old_local_;
  std::vector<double> precursor_new_local_;

  SetSourceFunction active_set_source_function_;

  bool initialized_ = false;

  /// Data carriers needed to run the sweep on GPU.
  std::shared_ptr<TotalXSCarrier> total_xs_carrier_ = nullptr;
  std::shared_ptr<OutflowCarrier> outflow_carrier_ = nullptr;
  std::shared_ptr<MeshCarrier> mesh_carrier_ = nullptr;

  /// Memory pinners for source moments and destination phi.
  std::shared_ptr<MemoryPinner<double>> source_pinner_ = nullptr;
  std::shared_ptr<MemoryPinner<double>> phi_pinner_ = nullptr;

  /// Flag indicating if GPU acceleration is enabled.
  bool use_gpus_;

private:
  void InitializeRuntimeCore();
  void ValidateRuntimeModeConfiguration() const;
  void InitializeSources();
  /// Initializes parallel arrays.
  void InitializeParrays();
  std::string MakeScalarFluxFieldFunctionName(unsigned int g, unsigned int m) const;
  std::vector<double> ComputeScalarFluxFieldFunctionData(unsigned int g, unsigned int m) const;
  double ComputeFieldFunctionPowerScaleFactor(double power_normalization_target) const;
  const std::vector<double>* GetFieldFunctionCoefficients(const MultiGroupXS& xs,
                                                          const std::string& xs_name) const;
  std::vector<double> ComputeXSFieldFunctionData(const std::string& xs_name) const;
  std::vector<double> ComputePowerFieldFunctionData(double& local_total_power) const;
  /// Initializes data carriers to GPUs and memory pinner.
  void InitializeGPUExtras();
  /// Reset data carriers to null and unpin memory.
  void ResetGPUCarriers();

  /// Initialize groupsets
  void InitializeGroupsets(const InputParameters& params);

  /// Initializes materials
  void InitializeXSMap(const InputParameters& params);
  void InitializeMaterials();

  /// Initialize sources
  void InitializeSources(const InputParameters& params);

  void ParseOptions(const InputParameters& input);

  static std::filesystem::path BuildRestartPath(const std::string& path_stem);

  /// Checks if the current CPU is associated with any GPU.
  static void CheckCapableDevices();

public:
  /// Max number of DOFs per cell that the sweep kernel on GPU can handle.
  static constexpr std::uint32_t max_dofs_gpu = 10;

  /// Returns the input parameters for this object.
  static InputParameters GetInputParameters();

  static InputParameters GetOptionsBlock();

  static InputParameters GetXSMapEntryBlock();
};

} // namespace opensn
