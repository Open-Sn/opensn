// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/point_source/point_source.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/volumetric_source/volumetric_source.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include <petscksp.h>
#include <any>
#include <chrono>

namespace opensn
{

class MPICommunicatorSet;
class GridFaceHistogram;
class TimeIntegration;
class AGSLinearSolver;
class WGSLinearSolver;
struct WGSContext;

/// Base class for all Linear Boltzmann Solvers.
class LBSProblem : public Problem
{
public:
  explicit LBSProblem(const std::string& name, std::shared_ptr<MeshContinuum> grid);

  /// Input parameters based construction.
  explicit LBSProblem(const InputParameters& params);

  LBSProblem(const LBSProblem&) = delete;

  LBSProblem& operator=(const LBSProblem&) = delete;

  ~LBSProblem() override;

  /// Returns a reference to the solver options.
  LBSOptions& GetOptions();

  /// Returns a constant reference to the solver options.
  const LBSOptions& GetOptions() const;

  static InputParameters GetOptionsBlock();

  static InputParameters GetBoundaryOptionsBlock();

  static InputParameters GetXSMapEntryBlock();

  void SetOptions(const InputParameters& input);

  void SetBoundaryOptions(const InputParameters& params);

  void SetAdjoint(bool adjoint);

  /// Returns the number of moments for the solver. This will only be non-zero after initialization.
  size_t GetNumMoments() const;

  /// Returns the number of groups for the solver. This will only be non-zero after initialization.
  size_t GetNumGroups() const;

  /// Returns the scattering order for the solver. This will only be non-zero after initialization.
  unsigned int GetScatteringOrder() const;

  /**
   * Returns the number of precursors for the solver. This will only be non-zero after
   * initialization.
   */
  size_t GetNumPrecursors() const;

  /**
   * Returns the maximum number of precursors defined on any material. This will only be non-zero
   * after initialization.
   */
  size_t GetMaxPrecursorsPerMaterial() const;

  /**
   * Adds a group to the list of groups. If group id < 0, the id will be logically derived from the
   * list size. If >= 0 the id will be set to the id specified.
   */
  void AddGroup(int id);

  const std::vector<LBSGroup>& GetGroups() const;

  /**
   * Adds a groupset to the list of groupsets. The groupset id will be logically derived from the
   * list size.
   */
  void AddGroupset();

  std::vector<LBSGroupset>& GetGroupsets();

  const std::vector<LBSGroupset>& GetGroupsets() const;

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
  void ClearBoundaries();

  size_t& GetLastRestartTime();

  /// Returns a reference to the map of material ids to XSs.
  const std::map<int, std::shared_ptr<MultiGroupXS>>& GetMatID2XSMap() const;

  /// Obtains a reference to the grid.
  std::shared_ptr<MeshContinuum> GetGrid() const;

  /// Get pointer to carriers.
  void* GetCarrier(std::uint32_t idx) { return carriers_.at(idx); }

  /// Get pointer to pinners.
  void* GetPinner(std::uint32_t idx) { return pinners_.at(idx); }

  /// Obtains a reference to the spatial discretization.
  const SpatialDiscretization& GetSpatialDiscretization() const;

  /// Returns read-only access to the unit cell matrices.
  const std::vector<UnitCellMatrices>& GetUnitCellMatrices() const;

  /// Returns read-only access to the unit ghost cell matrices.
  const std::map<uint64_t, UnitCellMatrices>& GetUnitGhostCellMatrices() const;

  /// Returns a reference to the list of local cell transport views.
  const std::vector<CellLBSView>& GetCellTransportViews() const;

  /// Read/Write access to the boundary preferences.
  std::map<uint64_t, BoundaryPreference>& GetBoundaryPreferences();

  /// Obtains a reference to the unknown manager for flux-moments.
  const UnknownManager& GetUnknownManager() const;

  /// Returns the local node count for the flux-moments data structures.
  size_t GetLocalNodeCount() const;

  /// Returns the global node count for the flux-moments data structures.
  size_t GetGlobalNodeCount() const;

  /// Read/write access to source moments vector.
  std::vector<double>& GetQMomentsLocal();

  /// Read access to source moments vector.
  const std::vector<double>& GetQMomentsLocal() const;

  /// Read/write access to exterior src moments vector.
  std::vector<double>& GetExtSrcMomentsLocal();

  /// Read access to exterior src moments vector.
  const std::vector<double>& GetExtSrcMomentsLocal() const;

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

  /// Read/write access to the cell-wise densities.
  std::vector<double>& GetDensitiesLocal();

  /// Read access to the cell-wise densities.
  const std::vector<double>& GetDensitiesLocal() const;

  SetSourceFunction GetActiveSetSourceFunction() const;

  std::shared_ptr<AGSLinearSolver> GetAGSSolver();

  std::vector<std::shared_ptr<LinearSolver>>& GetWGSSolvers();

  WGSContext& GetWGSContext(int groupset_id);

  /**
   * Gets the local and global number of iterative unknowns. This normally is only the flux moments,
   * however, the sweep based solvers might include delayed angular fluxes in this number.
   */
  virtual std::pair<size_t, size_t> GetNumPhiIterativeUnknowns();

  /// Gets the local handle of a flux-moment based field function.
  size_t MapPhiFieldFunction(size_t g, size_t m) const;

  /// Returns the power generation field function, if enabled.
  std::shared_ptr<FieldFunctionGridBased> GetPowerFieldFunction() const;

  void Initialize() override;

  /// Initializes default materials and physics materials.
  void InitializeMaterials();

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

  /// Copy relevant section of phi_old to the field functions.
  void UpdateFieldFunctions();

  /// Sets the internal phi vector to the value in the associated field function.
  void SetPhiFromFieldFunctions(PhiSTLOption which_phi,
                                const std::vector<size_t>& m_indices,
                                const std::vector<size_t>& g_indices);

  /**
   * A method for post-processing an adjoint solution.
   *
   * @note This does nothing for diffusion-based solvers.
   */
  virtual void ReorientAdjointSolution() {};

protected:
  /// Performs general input checks before initialization continues.
  virtual void PerformInputChecks();

  /// Prints header information of simulation.
  void PrintSimHeader();

  virtual void InitializeSpatialDiscretization();

  void ComputeUnitIntegrals();

  /// Initializes common groupset items.
  void InitializeGroupsets();

  /// Computes the number of moments for the given mesher types
  void ValidateAndComputeScatteringMoments();

  /// Initializes parallel arrays.
  virtual void InitializeParrays();

  void InitializeFieldFunctions();

  /// Initializes boundaries.
  virtual void InitializeBoundaries() {}

  virtual void InitializeSolverSchemes();

  virtual void InitializeWGSSolvers() {};

  /// Initializes data carriers to GPUs and memory pinner.
  void InitializeGPUExtras();

  /// Reset data carriers to null and unpin memory.
  void ResetGPUCarriers();

  virtual void ZeroSolutions() = 0;

  LBSOptions options_;
  size_t num_moments_ = 0;
  size_t num_groups_ = 0;
  unsigned int scattering_order_ = 0;
  size_t num_precursors_ = 0;
  size_t max_precursors_per_material_ = 0;

  std::vector<LBSGroup> groups_;
  std::vector<LBSGroupset> groupsets_;

  std::map<int, std::shared_ptr<MultiGroupXS>> block_id_to_xs_map_;

  std::vector<std::shared_ptr<PointSource>> point_sources_;
  std::vector<std::shared_ptr<VolumetricSource>> volumetric_sources_;

  std::shared_ptr<MeshContinuum> grid_;
  std::shared_ptr<SpatialDiscretization> discretization_ = nullptr;

  std::vector<CellFaceNodalMapping> grid_nodal_mappings_;
  std::shared_ptr<MPICommunicatorSet> grid_local_comm_set_ = nullptr;

  std::vector<UnitCellMatrices> unit_cell_matrices_;
  std::map<uint64_t, UnitCellMatrices> unit_ghost_cell_matrices_;
  std::vector<CellLBSView> cell_transport_views_;

  std::map<uint64_t, BoundaryPreference> boundary_preferences_;

  UnknownManager flux_moments_uk_man_;

  size_t max_cell_dof_count_ = 0;
  size_t min_cell_dof_count_ = 0;
  uint64_t local_node_count_ = 0;
  uint64_t global_node_count_ = 0;

  std::vector<double> q_moments_local_, ext_src_moments_local_;
  std::vector<double> phi_new_local_, phi_old_local_;
  std::vector<double> precursor_new_local_;
  std::vector<double> densities_local_;

  SetSourceFunction active_set_source_function_;

  std::shared_ptr<AGSLinearSolver> ags_solver_;
  std::vector<std::shared_ptr<LinearSolver>> wgs_solvers_;

  std::map<std::pair<size_t, size_t>, size_t> phi_field_functions_local_map_;
  size_t power_gen_fieldfunc_local_handle_ = 0;

  /// Time integration parameter meant to be set by an executor
  std::shared_ptr<const TimeIntegration> time_integration_ = nullptr;

  /**
   * @brief Data carriers for necessary data to run the sweep on GPU.
   * @details These objects manage GPU memory allocation automatically, organize cross-section,
   * outflow, and mesh data into contiguous memory on the CPU, and handle copying it to the GPU.
   *
   * There are 3 carriers, respectively for cross sections, outflow and mesh.
   */
  std::array<void*, 3> carriers_ = {nullptr, nullptr, nullptr};

  /// Memory pinner for source moments and destination phi.
  std::array<void*, 2> pinners_ = {nullptr, nullptr};

  /// Flag indicating if GPU acceleration is enabled.
  bool use_gpus_;

  /// Checks if the current CPU is associated with any GPU.
  static void CheckCapableDevices();

public:
  static std::map<std::string, uint64_t> supported_boundary_names;
  static std::map<uint64_t, std::string> supported_boundary_ids;

  /// Returns the input parameters for this object.
  static InputParameters GetInputParameters();
};

} // namespace opensn
