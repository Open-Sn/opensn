// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn
{

/**
 * \brief CBC sweep chunk for device
 */
class CBCD_SweepChunk : public CBCSweepChunk
{
public:
  CBCD_SweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  /// Check to save angular fluxes.
  bool IsSavingAngularFluxes() const { return save_angular_flux_; }

  /// Get the groupset's group index.
  unsigned int GetGroupsetGroupIndex() const { return groupset_.groups.front().id; }

  /// Get cell transport view for a local cell.
  const CellLBSView& GetCellTransportView(std::uint64_t cell_local_id) const
  {
    return cell_transport_views_[cell_local_id];
  }

  /// Device-analogue of Sweep.
  void GPUSweep(AngleSet& angle_set, const std::vector<std::uint64_t>& cell_local_ids);

  /// Create vector of pointers to CBCD_AngleSet, CBCD_FLUDS, and CBCD_AngleSet crb::Streams
  /// in groupset's angle aggregation.
  void PopulateAngleSetsAndFLUDSLists();

  /// Get vector of angle sets.
  const std::vector<CBCD_AngleSet*>& GetAngleSets() const { return angle_sets_; }

  /// Get vector of FLUDS.
  const std::vector<CBCD_FLUDS*>& GetFLUDSList() const { return fluds_list_; }

  /// Get vector of streams.
  const std::vector<crb::Stream*>& GetStreamsList() const { return streams_list_; }

  /// Allocate storage for saved angular fluxes on device.
  void AllocateSavedPsiStorage();

  /// Copy phi and source moments to device.
  void CopyPhiAndSrcToDevice();

  /// Copy back outflow and phi from device.
  void CopyOutflowAndPhiBackToHost();

  /// Asynchronously copy saved psi from device.
  void CopySavedPsiFromDevice(AngleSet& angle_set);

  /// Copy saved psi to host.
  void CopySavedPsiBackToHost(AngleSet& angle_set);

private:
  DiscreteOrdinatesProblem& problem_;
  std::vector<CBCD_AngleSet*> angle_sets_;
  std::vector<CBCD_FLUDS*> fluds_list_;
  std::vector<crb::Stream*> streams_list_;
};

} // namespace opensn