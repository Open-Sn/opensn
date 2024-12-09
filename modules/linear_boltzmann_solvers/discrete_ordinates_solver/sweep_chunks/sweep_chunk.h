// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_structs.h"
#include <functional>

namespace opensn
{

/// Sweep work function
class SweepChunk
{
public:
  /// Convenient typdef for the moment call back function. See moment_callbacks.
  using MomentCallbackFunc = std::function<void(SweepChunk* sweeper, AngleSet* angle_set)>;

  /**
   * Functions of type MomentCallbackFunc can be added to the moment_callbacks
   * vector and these can be called from within functions taking a
   * LBSGroupset instance. The intention is that this function can
   * be used as a general interface to retrieve angular flux values
   */
  std::vector<MomentCallbackFunc> moment_callbacks;

  SweepChunk(std::vector<double>& destination_phi,
             std::vector<double>& destination_psi,
             const MeshContinuum& grid,
             const SpatialDiscretization& discretization,
             const std::vector<UnitCellMatrices>& unit_cell_matrices,
             std::vector<CellLBSView>& cell_transport_views,
             const std::vector<double>& densities,
             const std::vector<double>& source_moments,
             const LBSGroupset& groupset,
             const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
             int num_moments,
             int max_num_cell_dofs)
    : grid_(grid),
      discretization_(discretization),
      unit_cell_matrices_(unit_cell_matrices),
      cell_transport_views_(cell_transport_views),
      densities_(densities),
      source_moments_(source_moments),
      groupset_(groupset),
      xs_(xs),
      num_moments_(num_moments),
      max_num_cell_dofs_(max_num_cell_dofs),
      save_angular_flux_(not destination_psi.empty()),
      groupset_angle_group_stride_(groupset_.psi_uk_man_.GetNumberOfUnknowns() *
                                   groupset_.groups.size()),
      groupset_group_stride_(groupset_.groups.size()),
      destination_phi_(&destination_phi),
      destination_psi_(&destination_psi)
  {
  }

  /// Sweep chunks should override this.
  virtual void Sweep(AngleSet& angle_set) {}

  /// Sets the currently active FLUx Data Structure
  virtual void SetAngleSet(AngleSet& angle_set) {}

  /// For cell-by-cell methods or computing the residual on a single cell.
  virtual void SetCell(Cell const* cell_ptr, AngleSet& angle_set) {}

  virtual ~SweepChunk() = default;

protected:
  friend class SweepScheduler;

  /// Sets the location where flux moments are to be written.
  void SetDestinationPhi(std::vector<double>& phi) { destination_phi_ = (&phi); }

  /**
   * Sets the portion of the output flux moments vector corresponding to the groupset this
   * sweep chunk operates on to zero.
   */
  void ZeroDestinationPhi();

  /// Returns a reference to the output flux moments vector.
  std::vector<double>& GetDestinationPhi() { return *destination_phi_; }

  /// Sets the location where angular fluxes are to be written.
  void SetDestinationPsi(std::vector<double>& psi) { destination_psi_ = (&psi); }

  /// Sets all elements of the output angular flux vector to zero.
  void ZeroDestinationPsi() { (*destination_psi_).assign((*destination_psi_).size(), 0.0); }

  /// Returns a reference to the output angular flux vector.
  std::vector<double>& GetDestinationPsi() { return *destination_psi_; }

  /// Activates or deactives the surface src flag.
  void SetBoundarySourceActiveFlag(bool flag_value) { surface_source_active_ = flag_value; }

  /// Returns the surface src-active flag.
  bool IsSurfaceSourceActive() const { return surface_source_active_; }

  const MeshContinuum& grid_;
  const SpatialDiscretization& discretization_;
  const std::vector<UnitCellMatrices>& unit_cell_matrices_;
  std::vector<CellLBSView>& cell_transport_views_;
  const std::vector<double>& densities_;
  const std::vector<double>& source_moments_;
  const LBSGroupset& groupset_;
  const std::map<int, std::shared_ptr<MultiGroupXS>>& xs_;
  const int num_moments_;
  const int max_num_cell_dofs_;
  const bool save_angular_flux_;
  const size_t groupset_angle_group_stride_;
  const size_t groupset_group_stride_;

private:
  std::vector<double>* destination_phi_;
  std::vector<double>* destination_psi_;
  bool surface_source_active_ = false;
};

} // namespace opensn
