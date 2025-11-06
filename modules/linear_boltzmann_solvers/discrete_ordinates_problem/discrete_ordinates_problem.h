// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include <memory>

namespace opensn
{

/**
 * Base class for Discrete Ordinates solvers. This class mostly establishes utilities related to
 * sweeping. From here we can derive a steady-state, transient, adjoint, and k-eigenvalue solver.
 */
class DiscreteOrdinatesProblem : public LBSProblem
{
protected:
  using SweepOrderGroupingInfo = std::pair<UniqueSOGroupings, DirIDToSOMap>;

public:
  /// Static registration based constructor.
  explicit DiscreteOrdinatesProblem(const InputParameters& params);
  ~DiscreteOrdinatesProblem() override;

  const std::string& GetSweepType() const { return sweep_type_; }

  std::pair<size_t, size_t> GetNumPhiIterativeUnknowns() override;

  /// Read/write access to newest updated angular flux vector.
  std::vector<std::vector<double>>& GetPsiNewLocal();

  /// Read access to newest updated angular flux vector.
  const std::vector<std::vector<double>>& GetPsiNewLocal() const;

  void PrintSimHeader() override;

  void Initialize() override;

  /// Returns the sweep boundaries as a read only reference
  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& GetSweepBoundaries() const;

  /// Reorient an adjoint solution to account for backwards streaming.
  void ReorientAdjointSolution() override;

  /// Zeroes all the outflow data-structures required to compute balance.
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);

protected:
  explicit DiscreteOrdinatesProblem(const std::string& name,
                                    std::shared_ptr<MeshContinuum> grid_ptr);

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

  void ZeroSolutions() override;

  std::map<std::shared_ptr<AngularQuadrature>, SweepOrderGroupingInfo>
    quadrature_unq_so_grouping_map_;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::shared_ptr<SPDS>>>
    quadrature_spds_map_;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::unique_ptr<FLUDSCommonData>>>
    quadrature_fluds_commondata_map_;

  std::vector<int> verbose_sweep_angles_;
  const std::string sweep_type_;
  std::map<uint64_t, std::shared_ptr<SweepBoundary>> sweep_boundaries_;

  /// Max level size.
  std::size_t max_level_size_ = 0;
  /// Max angle-set size.
  std::size_t max_angleset_size_ = 0;
  /// Max group-set size.
  std::size_t max_groupset_size_ = 0;

  std::shared_ptr<GridFaceHistogram> grid_face_histogram_ = nullptr;

  std::vector<std::vector<double>> psi_new_local_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<DiscreteOrdinatesProblem> Create(const ParameterBlock& params);

private:
  void CreateFLUDSCommonDataForDevice();
  std::shared_ptr<FLUDS> CreateFLUDSForDevice(std::size_t num_groups,
                                              std::size_t num_angles,
                                              const FLUDSCommonData& common_data);

  /// Computes the number of moments for the given mesher types
  void ValidateAndComputeScatteringMoments();

  /**
   * This routine groups angle-indices to groups sharing the same sweep ordering. It also takes
   * geometry into account.
   */
  std::pair<UniqueSOGroupings, DirIDToSOMap>
  AssociateSOsAndDirections(std::shared_ptr<MeshContinuum> grid,
                            const AngularQuadrature& quadrature,
                            AngleAggregationType agg_type,
                            GeometryType lbs_geo_type);
};

} // namespace opensn
