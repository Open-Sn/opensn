// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/outflow/cell_outflow_view.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/flux_moment_layout.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_view.h"
#include <algorithm>
#include <functional>
#include <vector>

namespace opensn
{

/// Sweep work function
class SweepChunk
{
public:
  /// Convenient typdef for the moment call back function. See moment_callbacks.
  using MomentCallbackFunc = std::function<void(SweepChunk* sweeper, AngleSet* angle_set)>;

  SweepChunk(std::vector<double>& destination_phi,
             std::vector<double>& destination_psi,
             const std::shared_ptr<Mesh>& grid,
             const SpatialDiscretization& discretization,
             const std::vector<UnitCellMatrices>& unit_cell_matrices,
             const std::vector<CellLBSView>& cell_transport_views,
             std::vector<CellOutflowView>& cell_outflow_views,
             const std::vector<double>& source_moments,
             const LBSGroupset& groupset,
             const BlockID2XSMap& xs,
             unsigned int num_moments,
             unsigned int max_num_cell_dofs,
             unsigned int min_num_cell_dofs)
    : grid_(grid),
      discretization_(discretization),
      unit_cell_matrices_(unit_cell_matrices),
      cell_transport_views_(cell_transport_views),
      cell_outflow_views_(cell_outflow_views),
      source_moments_(source_moments),
      groupset_(groupset),
      xs_(xs),
      num_moments_(num_moments),
      max_num_cell_dofs_(max_num_cell_dofs),
      min_num_cell_dofs_(min_num_cell_dofs),
      groupset_angle_group_stride_(groupset_.psi_uk_man_.GetNumberOfUnknowns() *
                                   groupset_.GetNumGroups()),
      groupset_group_stride_(groupset_.GetNumGroups()),
      destination_phi_(destination_phi),
      destination_psi_(destination_psi),
      cell_(nullptr)
  {
  }

  /// Sweep chunks should override this.
  virtual void Sweep(AngleSet& angle_set);

  /// Sets the currently active angleset.
  virtual void SetAngleSet(AngleSet& angle_set);

  /// Sets the active cell for cell-by-cell sweeps.
  void SetCell(Cell const* cell_ptr);

  /**
   * Zero the portion of the output flux moments vector corresponding to the groupset for this
   * sweep chunk.
   */
  virtual void ZeroDestinationPhi();

  /// Sets all elements of the output angular flux vector to zero.
  void ZeroDestinationPsi() { (destination_psi_).assign((destination_psi_).size(), 0.0); }

  /// Returns true if angular flux storage is allocated.
  bool SaveAngularFluxEnabled() const { return not destination_psi_.empty(); }

  /// Activates or deactivates the surface src flag.
  void SetBoundarySourceActiveFlag(bool flag_value) { surface_source_active_ = flag_value; }

  /// Returns the surface src-active flag.
  bool IsSurfaceSourceActive() const { return surface_source_active_; }

  /// Include contribution from previous psi_n in RHS source term
  void IncludeRHSTimeTerm(bool status) { include_rhs_time_term_ = status; }

  virtual bool IsTimeDependent() const { return false; }

  /**
   * Overwrite this groupset's scalar flux moments and outflows with sums over the stored angular
   * flux in a fixed order: angle sets in aggregation order, then their directions.
   *
   * \param angle_agg Angle aggregation of the groupset; its order fixes the summation order.
   */
  virtual void SetDestinationsFromAngularFlux(const AngleAggregation& angle_agg);

  /**
   * Prepares the chunk for a sweep.
   * The sweep scheduler calls this before the first angle set of a sweep, after the destination
   * flux moments of the groupset are zeroed. A chunk with compact groupset moments copies the
   * groupset's source moments here.
   */
  virtual void BeginSweep();

  /**
   * Finishes a sweep of the chunk.
   * The sweep scheduler calls this after the last angle set of a sweep. A chunk with compact
   * groupset moments copies the groupset's flux moments to the destination here.
   */
  virtual void EndSweep();

  /**
   * Where a sweep reads the groupset's source moments and adds its flux moments.
   * Moment m of group g of the groupset at node n is entry n * node_stride + m * moment_stride +
   * group_offset + g of both arrays. Initialize it with braces, which zero any members left out.
   */
  struct GroupsetMoments
  {
    const double* source;
    double* phi;
    std::size_t node_stride;
    std::size_t moment_stride;
    std::size_t group_offset;

    /// Index of moment m of the groupset's first group at node n.
    std::size_t Index(std::size_t n, unsigned int m) const
    {
      return n * node_stride + m * moment_stride + group_offset;
    }
  };

  virtual ~SweepChunk() = default;

  /**
   * Functions of type MomentCallbackFunc can be added to the moment_callbacks
   * vector and these can be called from within functions taking a
   * LBSGroupset instance. The intention is that this function can
   * be used as a general interface to retrieve angular flux values
   */
  std::vector<MomentCallbackFunc> moment_callbacks;

protected:
  /**
   * Makes the sweeps of this chunk use contiguous copies of the groupset's moments.
   * When `layout` does not store the groupset's moments of a node together, that is, when the
   * moment stride of the groupset exceeds its group count, a sweep reads the groupset's source
   * moments from, and adds its flux moments to, contiguous copies of the groupset's part of the
   * moments. The moments of a node then lie together, so a kernel reads them in fewer memory
   * accesses. The copies exist only during a sweep. A chunk enables this only if every kernel it
   * runs reads the moments through GetGroupsetMoments(). `layout` is the layout of the source and
   * destination moments.
   */
  void EnableCompactGroupsetMoments(const FluxMomentLayout& layout);

  /// Where the current sweep reads the source moments and adds the flux moments.
  GroupsetMoments GetGroupsetMoments();

  /**
   * Overwrite the scalar flux moments and outflows with fixed-order sums over an angular flux.
   *
   * \param angle_agg Angle aggregation of the groupset; its order fixes the summation order.
   * \param psi Angular flux laid out like destination_psi_.
   */
  void ComputeDestinationsFromAngularFlux(const AngleAggregation& angle_agg,
                                          const std::vector<double>& psi);

  /**
   * Overwrite the scalar flux moments and outflows of a time-dependent sweep.
   *
   * The moments are those of psi^{n+theta} = theta psi^{n+1} + (1 - theta) psi^n, while
   * destination_psi_ holds psi^{n+1}.
   *
   * \param angle_agg Angle aggregation of the groupset; its order fixes the summation order.
   * \param psi_old Angular flux psi^n at the start of the time step.
   * \param theta Time-integration parameter of the theta scheme.
   */
  void ComputeTimeDependentDestinationsFromAngularFlux(const AngleAggregation& angle_agg,
                                                       const std::vector<double>& psi_old,
                                                       double theta);

  const std::shared_ptr<Mesh> grid_;
  const SpatialDiscretization& discretization_;
  const std::vector<UnitCellMatrices>& unit_cell_matrices_;
  const std::vector<CellLBSView>& cell_transport_views_;
  std::vector<CellOutflowView>& cell_outflow_views_;
  const std::vector<double>& source_moments_;
  const LBSGroupset& groupset_;
  const BlockID2XSMap& xs_;
  const unsigned int num_moments_;
  const unsigned int max_num_cell_dofs_;
  const unsigned int min_num_cell_dofs_;
  const size_t groupset_angle_group_stride_;
  const size_t groupset_group_stride_;
  std::vector<double>& destination_phi_;
  std::vector<double>& destination_psi_;
  /// Angular flux at time n + theta, used by ComputeTimeDependentDestinationsFromAngularFlux.
  std::vector<double> psi_theta_;
  bool surface_source_active_ = false;
  bool include_rhs_time_term_ = true;
  const Cell* cell_;

private:
  /// Number of (node, moment) pairs of the moments; zero until EnableCompactGroupsetMoments.
  std::size_t num_node_moments_ = 0;
  /// Index of node 0, moment 0 of the groupset's first group in the moments.
  std::size_t groupset_base_ = 0;
  /// Distance between consecutive moments of one group, from the flux-moment layout.
  std::size_t moment_stride_ = 0;
  bool compact_groupset_moments_ = false;
  std::vector<double> compact_source_;
  std::vector<double> compact_phi_;
};

/**
 * Evaluate the angular source of one direction at the nodes of a cell for every group of a
 * groupset: nodal_source[i * gs_size + g] = sum over m of m2d_row[m] * q(i, m, gs_gi + g).
 *
 * The source moments of a node are contiguous in memory, ordered by moment and then by group, so
 * this reads them in memory order. The moments are summed in increasing order. A sweep kernel
 * calls this once per cell and direction and then applies the mass matrix once per group, not
 * once per moment. `nodal_source` must hold num_nodes * gs_size values.
 */
inline void
ComputeAngularSource(const CellLBSView& cell_transport_view,
                     std::size_t num_nodes,
                     unsigned int num_moments,
                     unsigned int gs_gi,
                     std::size_t gs_size,
                     const double* m2d_row,
                     const std::vector<double>& source_moments,
                     double* nodal_source)
{
  std::fill_n(nodal_source, num_nodes * gs_size, 0.0);
  for (std::size_t i = 0; i < num_nodes; ++i)
  {
    double* __restrict s_i = nodal_source + i * gs_size;
    for (unsigned int m = 0; m < num_moments; ++m)
    {
      const double w = m2d_row[m];
      const double* __restrict q = &source_moments[cell_transport_view.MapDOF(i, m, gs_gi)];
      for (std::size_t g = 0; g < gs_size; ++g)
        s_i[g] += w * q[g];
    }
  }
}

} // namespace opensn
