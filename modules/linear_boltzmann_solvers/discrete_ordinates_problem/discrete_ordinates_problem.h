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
    Default = 0,
    SteadyState = 1,
    TimeDependent = 2
  };

  void SetSweepChunkMode(SweepChunkMode mode);

  void ResetSweepChunkMode() { sweep_chunk_mode_.reset(); }

  bool IsTimeDependent() const override
  {
    return sweep_chunk_mode_.value_or(SweepChunkMode::Default) == SweepChunkMode::TimeDependent;
  }

  void SetTimeDependentMode() override;

  void SetSteadyStateMode() override;

  std::shared_ptr<SweepChunk> CreateSweepChunk(LBSGroupset& groupset)
  {
    return SetSweepChunk(groupset);
  }

  void ResetMode(SweepChunkMode target_mode);

  /// Rebuild WGS/AGS solver schemes (e.g., after changing sweep chunk mode).
  void ReinitializeSolverSchemes();

  ~DiscreteOrdinatesProblem() override;

  using BoundaryDefinition = std::pair<LBSBoundaryType, std::shared_ptr<SweepBoundary>>;

  const std::string& GetSweepType() const { return sweep_type_; }

  std::pair<std::uint64_t, std::uint64_t> GetNumPhiIterativeUnknowns() override;

  /// Read/write access to newest updated angular flux vector.
  std::vector<std::vector<double>>& GetPsiNewLocal();

  /// Read access to newest updated angular flux vector.
  const std::vector<std::vector<double>>& GetPsiNewLocal() const;

  /// Read/write access to newest updated angular flux vector.
  std::vector<std::vector<double>>& GetPsiOldLocal();

  /// Read access to previous angular flux vector.
  const std::vector<std::vector<double>>& GetPsiOldLocal() const;

  size_t GetMaxLevelSize() const;

  size_t GetMaxGroupsetSize() const;

  size_t GetMaxAngleSetSize() const;

  /// Copy psi_new to psi_old
  void UpdatePsiOld() override;

  void PrintSimHeader() override;

  /// Returns the sweep boundaries as a read only reference
  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& GetSweepBoundaries() const;

  const std::map<uint64_t, BoundaryDefinition>& GetBoundaryDefinitions() const;

  /// Reorient an adjoint solution to account for backwards streaming.
  void ReorientAdjointSolution() override;

  /// Zeroes all the outflow data-structures required to compute balance.
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);

  /**
   * Create (if needed) and return angular flux field functions for the given groups and angles.
   *
   * Angles are indices into the groupset quadrature associated with each group.
   */
  std::vector<std::shared_ptr<FieldFunctionGridBased>>
  GetAngularFluxFieldFunctionList(const std::vector<unsigned int>& groups,
                                  const std::vector<size_t>& angles);

  void SetSaveAngularFlux(bool save) override;

  /// Update angular flux field functions from psi_new_local_.
  void UpdateAngularFieldFunctions();

  void SetBoundaryOptions(const InputParameters& params) override;
  void ClearBoundaries() override;

  void CopyPhiAndSrcToDevice();
  void CopyPhiAndOutflowBackToHost();

protected:
  /// Factory-only constructor.
  explicit DiscreteOrdinatesProblem(const InputParameters& params);

  /// Build sweep/runtime data structures once base runtime is available.
  void BuildRuntime();

  void InitializeBoundaries() override;

  /// Initializes Within-GroupSet solvers.
  void InitializeWGSSolvers() override;

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

  /// Clears all the sweep orderings for a groupset in preperation for another.
  void ResetSweepOrderings(LBSGroupset& groupset);

  /// Sets up the sweek chunk for the given discretization method.
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);

  void ZeroPsi() override;

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
  bool forced_save_angular_flux_for_transient_ = false;
  bool save_angular_flux_before_transient_ = false;

  std::map<std::tuple<size_t, size_t, size_t>, size_t> angular_flux_field_functions_local_map_;

private:
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

public:
  static InputParameters GetInputParameters();
  static InputParameters GetBoundaryOptionsBlock();
  static std::shared_ptr<DiscreteOrdinatesProblem> Create(const ParameterBlock& params);
};

} // namespace opensn
