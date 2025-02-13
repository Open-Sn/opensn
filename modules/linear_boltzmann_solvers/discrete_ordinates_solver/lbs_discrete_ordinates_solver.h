// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"
#include <memory>

namespace opensn
{

/**
 * Base class for Discrete Ordinates solvers. This class mostly establishes utilities related to
 * sweeping. From here we can derive a steady-state, transient, adjoint, and k-eigenvalue solver.
 */
class DiscreteOrdinatesSolver : public LBSSolver
{
protected:
  using SweepOrderGroupingInfo = std::pair<UniqueSOGroupings, DirIDToSOMap>;

public:
  /// Static registration based constructor.
  explicit DiscreteOrdinatesSolver(const InputParameters& params);
  ~DiscreteOrdinatesSolver() override;

  const std::string& GetSweepType() const { return sweep_type_; }

  std::pair<size_t, size_t> GetNumPhiIterativeUnknowns() override;
  void Initialize() override;

  /// Reorient an adjoint solution to account for backwards streaming.
  void ReorientAdjointSolution() override;

  /// Zeroes all the outflow data-structures required to compute balance.
  void ZeroOutflowBalanceVars(LBSGroupset& groupset);

  /// Compute balance
  void ComputeBalance();

  /**
   * Computes the angular flux based leakage from boundary surfaces.
   * \param groupset_id The groupset for which to compute the leakage.
   * \param boundary_id The boundary id for which to perform the integration.
   *
   * \return The leakage as a value.
   */
  std::vector<double> ComputeLeakage(unsigned int groupset_id, uint64_t boundary_id) const;

  /**
   * Computes the group-wise angular flux-based leakage from the specified boundaries.
   *
   * \param boundary_ids The boundary ids to compute leakages on.
   * \return A map of boundary ids to group-wise leakages.
   */
  std::map<uint64_t, std::vector<double>>
  ComputeLeakage(const std::vector<uint64_t>& boundary_ids) const;

protected:
  explicit DiscreteOrdinatesSolver(const std::string& name,
                                   std::shared_ptr<MeshContinuum> grid_ptr);

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
   *
   * The Template FLUDS can be scaled with number of angles and groups which
   * provides us with the angle-set-subset- and groupset-subset capability.
   */
  void InitializeSweepDataStructures();

  /// Initializes fluds_ data structures.
  void InitFluxDataStructures(LBSGroupset& groupset);

  /// Clears all the sweep orderings for a groupset in preperation for another.
  void ResetSweepOrderings(LBSGroupset& groupset);

  /// Sets up the sweek chunk for the given discretization method.
  virtual std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset);

  std::map<std::shared_ptr<AngularQuadrature>, SweepOrderGroupingInfo>
    quadrature_unq_so_grouping_map_;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::shared_ptr<SPDS>>>
    quadrature_spds_map_;
  std::map<std::shared_ptr<AngularQuadrature>, std::vector<std::unique_ptr<FLUDSCommonData>>>
    quadrature_fluds_commondata_map_;

  std::vector<size_t> verbose_sweep_angles_;
  const std::string sweep_type_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<DiscreteOrdinatesSolver> Create(const ParameterBlock& params);

protected:
  /**
   * This routine groups angle-indices to groups sharing the same sweep ordering. It also takes
   * geometry into account.
   */
  static std::pair<UniqueSOGroupings, DirIDToSOMap>
  AssociateSOsAndDirections(const MeshContinuum& grid,
                            const AngularQuadrature& quadrature,
                            AngleAggregationType agg_type,
                            GeometryType lbs_geo_type);
};

} // namespace opensn
