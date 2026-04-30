// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/parameters/parameter_block.h"
#include <memory>
#include <optional>
#include <tuple>

namespace opensn
{
class FieldFunctionGridBased;
class AGSLinearSolver;
class LinearSolver;
struct BalanceTable;
struct WGSContext;

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

  /// Construction and transport-mode controls.
  bool IsTimeDependent() const override
  {
    return sweep_chunk_mode_.value_or(SweepChunkMode::DEFAULT) == SweepChunkMode::TIME_DEPENDENT;
  }

  void SetTimeDependentMode() override;

  void SetSteadyStateMode() override;

  ~DiscreteOrdinatesProblem() override;

  using BoundaryDefinition = std::pair<LBSBoundaryType, std::shared_ptr<SweepBoundary>>;

  /// Problem metadata and solver access.
  const std::string& GetSweepType() const { return sweep_type_; }

  std::pair<size_t, size_t> GetNumPhiIterativeUnknowns() override;

  std::shared_ptr<AGSLinearSolver> GetAGSSolver();

  std::shared_ptr<LinearSolver> GetWGSSolver(size_t groupset_id);
  size_t GetNumWGSSolvers();

  WGSContext& GetWGSContext(int groupset_id);

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

  /// Balance and output helpers.
  BalanceTable ComputeBalanceTable(double scaling_factor = 1.0);

  void ComputeBalance(double scaling_factor = 1.0);

  void PrintSimHeader() override;

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

  /**
   * Supported runtime discrete-ordinates reconfiguration.
   *
   * These methods update the dependent sweep, boundary, acceleration, and solver data
   * owned by the problem. Lower-level sweep chunk and solver-scheme controls remain
   * protected implementation details.
   */
  void SetSaveAngularFlux(bool save);

  void SetBlockID2XSMap(const BlockID2XSMap& xs_map) override;

  void SetBoundaryOptions(const InputParameters& params) override;
  void ClearBoundaries() override;

  void CopyPhiAndSrcToDevice();
  void CopyPhiAndOutflowBackToHost();

protected:
  /// Factory-only constructor.
  explicit DiscreteOrdinatesProblem(const InputParameters& params);

  /// Internal factory step: build sweep/runtime data once base runtime data is available.
  void BuildRuntime();

  void InitializeBoundaries() override;

  void InitializeSolverSchemes();
  /// Rebuild WGS/AGS solver schemes after runtime configuration changes.
  void ReinitializeSolverSchemes();

  void SetSweepChunkMode(SweepChunkMode mode);
  void ResetSweepChunkMode() { sweep_chunk_mode_.reset(); }
  void ResetMode(SweepChunkMode target_mode);

  void InitializeWGSContexts();

  /// Initializes Within-GroupSet solvers.
  void InitializeWGSSolvers();

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

  std::shared_ptr<SweepChunk> CreateSweepChunk(LBSGroupset& groupset)
  {
    return SetSweepChunk(groupset);
  }

  bool ReadProblemRestartData(hid_t file_id) override;
  bool WriteProblemRestartData(hid_t file_id) const override;
  void ResetDerivedSolutionVectors() override;
  void RebuildBoundaryRuntimeData();

  BoundaryDefinition CreateBoundaryFromParams(const InputParameters& params) const;
  std::shared_ptr<SweepBoundary> CreateSweepBoundary(uint64_t boundary_id) const;

  std::map<std::shared_ptr<AngularQuadrature>, SweepOrderGroupingInfo>
    quadrature_unq_so_grouping_map_;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::shared_ptr<SPDS>>>
    quadrature_spds_map_;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::unique_ptr<FLUDSCommonData>>>
    quadrature_fluds_commondata_map_;

  std::vector<int> verbose_sweep_angles_;
  const std::string sweep_type_;
  std::map<uint64_t, std::shared_ptr<SweepBoundary>> sweep_boundaries_;
  std::map<uint64_t, BoundaryDefinition> boundary_definitions_;
  std::optional<ParameterBlock> boundary_conditions_block_;

  /// Max level size.
  std::size_t max_level_size_ = 0;
  /// Max angle-set size.
  std::size_t max_angleset_size_ = 0;
  /// Max group-set size.
  unsigned int max_groupset_size_ = 0;

  std::shared_ptr<GridFaceHistogram> grid_face_histogram_ = nullptr;

  std::vector<std::vector<double>> psi_new_local_;
  std::vector<std::vector<double>> psi_old_local_;
  std::optional<SweepChunkMode> sweep_chunk_mode_;
  std::shared_ptr<AGSLinearSolver> ags_solver_;
  std::vector<std::shared_ptr<WGSContext>> wgs_contexts_;
  std::vector<std::shared_ptr<LinearSolver>> wgs_solvers_;

private:
  std::string
  MakeAngularFieldFunctionName(size_t groupset_id, unsigned int group, size_t angle) const;
  std::vector<double>
  ComputeAngularFieldFunctionData(size_t groupset_id, unsigned int group, size_t angle) const;

  void CreateAAHD_FLUDSCommonData();
  std::shared_ptr<FLUDS> CreateAAHD_FLUDS(unsigned int num_groups,
                                          std::size_t num_angles,
                                          const FLUDSCommonData& common_data);
  std::shared_ptr<AngleSet>
  CreateAAHD_AngleSet(size_t id,
                      unsigned int num_groups,
                      const SPDS& spds,
                      std::shared_ptr<FLUDS>& fluds,
                      std::vector<size_t>& angle_indices,
                      std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                      int maximum_message_size,
                      const MPICommunicatorSet& in_comm_set);
  std::shared_ptr<SweepChunk> CreateAAHD_SweepChunk(LBSGroupset& groupset);

  void CreateCBCD_FLUDSCommonData();
  std::shared_ptr<FLUDS> CreateCBCD_FLUDS(std::size_t num_groups,
                                          std::size_t num_angles,
                                          std::size_t num_local_cells,
                                          const FLUDSCommonData& common_data,
                                          const UnknownManager& psi_uk_man,
                                          const SpatialDiscretization& sdm,
                                          bool save_angular_flux);

  std::shared_ptr<AngleSet>
  CreateCBCD_AngleSet(size_t id,
                      size_t num_groups,
                      const SPDS& spds,
                      std::shared_ptr<FLUDS>& fluds,
                      std::vector<size_t>& angle_indices,
                      std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                      const MPICommunicatorSet& in_comm_set);
  std::shared_ptr<SweepChunk> CreateCBCDSweepChunk(LBSGroupset& groupset);

  /**
   * This routine groups angle-indices to groups sharing the same sweep ordering. It also takes
   * geometry into account.
   */
  std::pair<UniqueSOGroupings, DirIDToSOMap>
  AssociateSOsAndDirections(std::shared_ptr<MeshContinuum> grid,
                            const AngularQuadrature& quadrature,
                            AngleAggregationType agg_type,
                            GeometryType lbs_geo_type);

  void UpdateAngularFluxStorage();
  void UpdateAngularFluxFieldFunction(FieldFunctionGridBased& ff,
                                      size_t groupset_id,
                                      unsigned int group,
                                      size_t angle);

public:
  static InputParameters GetInputParameters();
  static InputParameters GetBoundaryOptionsBlock();
  static std::shared_ptr<DiscreteOrdinatesProblem> Create(const ParameterBlock& params);
};

} // namespace opensn
