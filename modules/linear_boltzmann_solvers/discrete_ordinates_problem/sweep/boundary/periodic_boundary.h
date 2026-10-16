// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/mesh/mesh_continuum/cell.h"
#include <cstdint>
#include <map>
#include <memory>
#include <vector>

namespace opensn
{

class MeshContinuum;

/// A Cartesian periodic face. Opposite-face flux is a delayed sweep unknown.
class PeriodicBoundary : public SweepBoundary
{
public:
  PeriodicBoundary(BoundaryBank& bank,
                   std::uint64_t boundary_id,
                   const std::shared_ptr<MeshContinuum>& grid,
                   const std::vector<LBSGroupset>& groupsets);

  void SetPartner(PeriodicBoundary& partner) { partner_ = &partner; }
  void GatherNodeGeometry();
  void MatchPartnerNodes();

  bool HasDelayedAngularFlux() const override { return true; }
  void InitializeAngleDependent(const std::vector<LBSGroupset>& groupsets) override;
  void ZeroOpposingDelayedAngularFluxOld(int groupset_id) override;
  size_t CountDelayedAngularDOFsNew(int groupset_id) const override;
  size_t CountDelayedAngularDOFsOld(int groupset_id) const override;
  void PrepareNewDelayedAngularFlux(int groupset_id) override;
  void AppendNewDelayedAngularDOFsToVector(int groupset_id,
                                           std::vector<double>& output) const override;
  void AppendOldDelayedAngularDOFsToVector(int groupset_id,
                                           std::vector<double>& output) const override;
  void AppendNewDelayedAngularDOFsToArray(int groupset_id,
                                          int64_t& index,
                                          double* buffer) const override;
  void AppendOldDelayedAngularDOFsToArray(int groupset_id,
                                          int64_t& index,
                                          double* buffer) const override;
  void
  SetNewDelayedAngularDOFsFromArray(int groupset_id, int64_t& index, const double* buffer) override;
  void
  SetOldDelayedAngularDOFsFromArray(int groupset_id, int64_t& index, const double* buffer) override;
  void SetNewDelayedAngularDOFsFromVector(int groupset_id,
                                          const std::vector<double>& values,
                                          size_t& index) override;
  void SetOldDelayedAngularDOFsFromVector(int groupset_id,
                                          const std::vector<double>& values,
                                          size_t& index) override;
  void CopyDelayedAngularFluxOldToNew(int groupset_id) override;
  void CopyDelayedAngularFluxNewToOld(int groupset_id) override;

  double* PsiIncoming(std::uint32_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int groupset_id,
                      unsigned int group_idx) override;
  double* PsiOutgoing(uint64_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int groupset_id) override;
  std::uint64_t
  GetOffsetToAngleset(const FaceNode& face_node, AngleSet& angleset, bool is_outgoing) override;

private:
  struct GroupsetData
  {
    std::vector<std::uint64_t> incoming_angle_index;
    std::uint64_t num_angles = 0;
    std::uint64_t num_incoming_angles = 0;
    std::uint64_t incoming_stride = 0;
    std::uint64_t outgoing_stride = 0;
  };

  std::vector<double> GatherOutgoing(int groupset_id) const;
  void AssembleIncoming(int groupset_id);
  double* IncomingNew(int groupset_id);
  double* IncomingOld(int groupset_id);
  const double* IncomingNew(int groupset_id) const;
  const double* IncomingOld(int groupset_id) const;
  size_t IncomingSize(int groupset_id) const;

  std::uint64_t boundary_id_;
  PeriodicBoundary* partner_ = nullptr;
  std::map<FaceNode, std::uint64_t> node_to_index_;
  // Each record is [node xyz, face-centroid xyz], in local face-node order.
  std::vector<double> local_geometry_;
  std::vector<double> global_geometry_;
  std::vector<std::uint64_t> global_node_offsets_;
  // Local incoming node -> global outgoing node on the partner boundary.
  std::vector<std::uint64_t> partner_node_index_;
  std::vector<GroupsetData> groupset_data_;
};

} // namespace opensn
