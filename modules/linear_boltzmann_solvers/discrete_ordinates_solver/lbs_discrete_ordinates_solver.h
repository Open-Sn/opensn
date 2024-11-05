// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

namespace opensn
{

class CBC_ASynchronousCommunicator;

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

  const std::string& SweepType() const { return sweep_type_; }

  std::pair<size_t, size_t> NumPhiIterativeUnknowns() override;
  void Initialize() override;
  void ScalePhiVector(PhiSTLOption which_phi, double value) override;
  void SetGSPETScVecFromPrimarySTLvector(const LBSGroupset& groupset,
                                         Vec x,
                                         PhiSTLOption which_phi) override;
  void SetPrimarySTLvectorFromGSPETScVec(const LBSGroupset& groupset,
                                         Vec x,
                                         PhiSTLOption which_phi) override;
  void GSScopedCopyPrimarySTLvectors(const LBSGroupset& groupset,
                                     PhiSTLOption from_which_phi,
                                     PhiSTLOption to_which_phi) override;
  void SetMultiGSPETScVecFromPrimarySTLvector(const std::vector<int>& groupset_ids,
                                              Vec x,
                                              PhiSTLOption which_phi) override;
  void SetPrimarySTLvectorFromMultiGSPETScVecFrom(const std::vector<int>& groupset_ids,
                                                  Vec x,
                                                  PhiSTLOption which_phi) override;
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
  explicit DiscreteOrdinatesSolver(const std::string& name);

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
