// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_definition.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_bank.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/parameters/parameter_block.h"
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace opensn
{

class FieldFunctionGridBased;
class AGSLinearSolver;
class LinearSolver;
struct BalanceTable;
struct WGSContext;
class AngularFluxFunction;
class BoundaryCarrier;

/**
 * Base class for discrete ordinates solvers.
 */
class DiscreteOrdinatesProblem : public LBSProblem
{
protected:
  using SweepOrderGroupingInfo = std::pair<UniqueSOGroupings, DirIDToSOMap>;

public:
  enum class SweepChunkMode
  {
    DEFAULT = 0,
    STEADY_STATE = 1,
    TIME_DEPENDENT = 2
  };

  /**
   * @name Construction and transport-mode controls
   * @{
   */
  bool IsTimeDependent() const override
  {
    return sweep_chunk_mode_.value_or(SweepChunkMode::DEFAULT) == SweepChunkMode::TIME_DEPENDENT;
  }

  void SetTimeDependentMode() override;

  void SetSteadyStateMode() override;

  void SetTime(double time) override;
  /** @} */

  ~DiscreteOrdinatesProblem() override;

  /**
   * @name Problem metadata and solver access
   * @{
   */
  const std::string& GetSweepType() const { return sweep_type_; }

  std::pair<std::uint64_t, std::uint64_t> GetNumPhiIterativeUnknowns() override;

  std::shared_ptr<AGSLinearSolver> GetAGSSolver();

  std::shared_ptr<LinearSolver> GetWGSSolver(size_t groupset_id);
  size_t GetNumWGSSolvers();

  WGSContext& GetWGSContext(int groupset_id);
  /** @} */

  /**
   * Internal angular-flux state access.
   *
   * These mutable references are used by sweep chunks, transient solvers, acceleration,
   * restart I/O, and angular-flux I/O. User-facing Python APIs expose copies or field
   * functions rather than mutable access to this storage.
   */
  /// Read/write access to newest updated angular flux vector.
  std::vector<std::vector<double>>& GetPsiNewLocal();

  /// Read access to newest updated angular flux vector.
  const std::vector<std::vector<double>>& GetPsiNewLocal() const;

  /// Read/write access to newest updated angular flux vector.
  std::vector<std::vector<double>>& GetPsiOldLocal();

  /// Read access to previous angular flux vector.
  const std::vector<std::vector<double>>& GetPsiOldLocal() const;

  void ZeroPsi();

  bool SaveAngularFluxEnabled() const { return options_.save_angular_flux; }

  size_t GetMaxLevelSize() const;

  size_t GetMaxGroupsetSize() const;

  size_t GetMaxAngleSetSize() const;

  /// Copy psi_new to psi_old
  void UpdatePsiOld();

  /**
   * @name Balance and output helpers
   * @{
   */
  BalanceTable ComputeBalanceTable(double scaling_factor = 1.0);

  void ComputeBalance(double scaling_factor = 1.0);

  /// Returns the sweep boundaries as a read-only reference.
  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& GetSweepBoundaries() const;

  const std::map<uint64_t, BoundaryDefinition>& GetBoundaryDefinitions() const;

  /// Reorient an adjoint solution to account for backwards streaming.
  void ReorientAdjointSolution() override;

  /// Zeroes all the outflow data-structures required to compute balance.
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);

  /**
   * Create angular flux field functions for the given groups and angles.
   *
   * Angles are indices into the groupset quadrature associated with each group.
   */
  std::vector<std::shared_ptr<FieldFunctionGridBased>>
  CreateAngularFluxFieldFunctionList(const std::vector<unsigned int>& groups,
                                     const std::vector<size_t>& angles);
  /** @} */

  /**
   * Supported runtime discrete-ordinates reconfiguration.
   *
   * These methods update the dependent sweep, boundary, acceleration, and solver data
   * owned by the problem. Lower-level sweep chunk and solver-scheme controls remain
   * protected implementation details.
   */
  void SetSaveAngularFlux(bool save);

  void SetBlockID2XSMap(const BlockID2XSMap& xs_map) override;

  void SetBoundaryOptions(const std::vector<InputParameters>& boundary_params, bool clear_existing);
  void ClearBoundaries() override;
  BoundaryCarrier* GetBoundaryCarrier() { return boundary_carrier_.get(); }
  const BoundaryCarrier* GetBoundaryCarrier() const { return boundary_carrier_.get(); }

  void CopyPhiAndSrcToDevice();
  void CopyPhiAndOutflowBackToHost();
  /**
   * Transfer data in boundary to device or vice-versa.
   * \param groupset_id Groupset ID.
   * \param host_to_device Flag indicating the direction of transfer.
   * \param force Force update.
   */
  void TransferDeviceBoundaryData(int groupset_id, bool host_to_device, bool force = false);

  bool HasUncollidedFlux() const { return not uncollided_flux_file_.empty(); }

  const std::vector<double>& GetFirstCollisionSourceMoments() const
  {
    return first_collision_source_moments_local_;
  }

  double GetUncollidedSourceRate() const { return uncollided_source_rate_; }
  double GetUncollidedOutflowRate() const { return uncollided_outflow_rate_; }

  /// Removes a previously applied uncollided scalar flux from the iterative state.
  /// Returns true when flux was removed.
  bool RemoveUncollidedFlux();

  /// Adds the uncollided scalar flux to the converged collided scalar flux.
  void ComputeFluxFromUncollided();

protected:
  /**
   * @name Construction and setup
   * @{
   */

  /// Factory-only constructor.
  explicit DiscreteOrdinatesProblem(const InputParameters& params);

  /// Internal factory step: build sweep/runtime data once base runtime data is available.
  void BuildRuntime();

  void InitializeBoundaries() override;
  /** @} */

  void PrintSimHeader() override;

  /**
   * @name Solver and sweep state
   * @{
   */

  void InitializeSolverSchemes();
  /// Rebuild WGS/AGS solver schemes after runtime configuration changes.
  void ReinitializeSolverSchemes();
  void ConfigureTransientSourceScopes();

  void SetSweepChunkMode(SweepChunkMode mode);
  void ResetSweepChunkMode() { sweep_chunk_mode_.reset(); }
  void ResetMode(SweepChunkMode target_mode);

  void InitializeBoundaryCarrier();

  /**
   * Sort the angle indices within each angle set so that one set maps exactly to another under all
   * reflected angle mappings.
   */
  void SortAngleSetsAngleIndices();

  void ResetBoundaryCarrier();

  /**
   * This routine initializes basic sweep datastructures that are agnostic of
   * the number of groups and essentially the groupsets. The routine rebuilds
   * the data structures i) `quadrature_unq_so_grouping_map_`,
   * ii) `quadrature_spds_map_` and iii) `quadrature_fluds_templates_map_`.
   * i) is a mapping, per quadrature, to a collection of angle-index-sets where
   * all the angles in a particular angleset share the same sweep ordering.
   * ii) is a mapping, per quadrature, to a collection of SPDSs where each
   * SPDS mirrors an angle-index-set in i)
   * iii) is again a mapping, per quadrature, to a collection of Template FLUDS
   * where each FLUDS mirrors a SPDS in ii).
   */
  void InitializeSweepDataStructures();

  /// Initializes fluds_ data structures.
  void InitFluxDataStructures(LBSGroupset& groupset);

  /// Clears all the sweep orderings for a groupset in preparation for another.
  void ResetSweepOrderings(LBSGroupset& groupset);

  /// Sets up the sweep chunk for the given discretization method.
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);
  /** @} */

  /**
   * @name Restart and runtime reconfiguration
   * @{
   */

  bool ReadProblemRestartData(hid_t file_id,
                              bool allow_transient_initialization_from_steady) override;
  bool WriteProblemRestartData(hid_t file_id) const override;
  void ResetDerivedSolutionVectors() override;
  void UpdateBoundaryDefinition(const InputParameters& params);
  void RebuildBoundaryRuntimeData();

  /**
   * @name Sweep dependency data
   * @{
   */
  std::map<std::shared_ptr<AngularQuadrature>, SweepOrderGroupingInfo>
    quadrature_unq_so_grouping_map_;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::shared_ptr<SPDS>>>
    quadrature_spds_map_;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::unique_ptr<FLUDSCommonData>>>
    quadrature_fluds_commondata_map_;

  const std::string sweep_type_;
  /** @} */

  /**
   * @name Boundary data
   * @{
   */
  BoundaryBank boundary_bank_;
  std::map<uint64_t, std::shared_ptr<SweepBoundary>> sweep_boundaries_;
  std::map<uint64_t, BoundaryDefinition> boundary_definitions_;
  std::shared_ptr<BoundaryCarrier> boundary_carrier_ = nullptr;
  std::optional<ParameterBlock> boundary_conditions_block_;
  bool boundary_runtime_data_initialized_ = false;
  bool has_reflecting_boundaries_ = false;
  bool has_time_dependent_boundaries_ = false;
  /** @} */

  /**
   * @name Sweep size metadata
   * @{
   */
  /// Max level size.
  std::size_t max_level_size_ = 0;
  /// Max angle-set size.
  std::size_t max_angleset_size_ = 0;
  /// Max group-set size.
  unsigned int max_groupset_size_ = 0;
  std::vector<std::vector<double>> psi_new_local_;
  std::vector<std::vector<double>> psi_old_local_;
  std::optional<SweepChunkMode> sweep_chunk_mode_;
  std::string uncollided_flux_file_;
  std::vector<double> uncollided_flux_moments_local_;
  std::vector<double> first_collision_source_moments_local_;
  double uncollided_source_rate_ = 0.0;
  double uncollided_outflow_rate_ = 0.0;
  bool uncollided_flux_applied_ = false;
  std::shared_ptr<AGSLinearSolver> ags_solver_;
  std::vector<std::shared_ptr<WGSContext>> wgs_contexts_;
  std::vector<std::shared_ptr<LinearSolver>> wgs_solvers_;
  /** @} */

  void InitializeFCS();

private:
  /**
   * @name Angular flux field-function helpers
   * @{
   */

  std::string
  MakeAngularFieldFunctionName(size_t groupset_id, unsigned int group, size_t angle) const;
  std::vector<double>
  ComputeAngularFieldFunctionData(size_t groupset_id, unsigned int group, size_t angle) const;

  void UpdateAAHD_FLUDSCommonDataWithBoundary();
  std::shared_ptr<FLUDS> CreateAAHD_FLUDS(unsigned int num_groups,
                                          std::size_t num_angles,
                                          const FLUDSCommonData& common_data);
  std::shared_ptr<AngleSet>
  CreateAAHD_AngleSet(size_t id,
                      const LBSGroupset& groupset,
                      const SPDS& spds,
                      std::shared_ptr<FLUDS>& fluds,
                      std::vector<size_t>& angle_indices,
                      std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                      int maximum_message_size,
                      const MPICommunicatorSet& in_comm_set);
  std::shared_ptr<SweepChunk> CreateAAHD_SweepChunk(LBSGroupset& groupset);

  std::shared_ptr<FLUDS> CreateCBCD_FLUDS(std::size_t num_groups,
                                          std::size_t num_angles,
                                          std::size_t num_local_cells,
                                          const FLUDSCommonData& common_data,
                                          const UnknownManager& psi_uk_man,
                                          const SpatialDiscretization& sdm,
                                          bool save_angular_flux);

  std::shared_ptr<AngleSet>
  CreateCBCD_AngleSet(size_t id,
                      const LBSGroupset& groupset,
                      const SPDS& spds,
                      std::shared_ptr<FLUDS>& fluds,
                      std::vector<size_t>& angle_indices,
                      std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                      const MPICommunicatorSet& in_comm_set);
  std::shared_ptr<SweepChunk> CreateCBCDSweepChunk(LBSGroupset& groupset);

  void UpdateAngularFluxStorage();
  void UpdateAngularFluxFieldFunction(FieldFunctionGridBased& ff,
                                      size_t groupset_id,
                                      unsigned int group,
                                      size_t angle);
  /** @} */

public:
  /**
   * @name Factory interface
   * @{
   */

  static InputParameters GetInputParameters();
  static InputParameters GetBoundaryOptionsBlock();
  static std::shared_ptr<DiscreteOrdinatesProblem> Create(const ParameterBlock& params);
  /** @} */
};

} // namespace opensn
