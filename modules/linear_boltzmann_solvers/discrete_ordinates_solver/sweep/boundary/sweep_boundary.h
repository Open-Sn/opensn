// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_structs.h"
#include <vector>

namespace opensn
{

/// Base class for sweep related boundaries.
class SweepBoundary
{
private:
  const LBSBoundaryType type_;
  const CoordinateSystemType coord_type_;
  /// Time value passed to boundary functions
  double evaluation_time_ = 0.0;

protected:
  std::vector<double> zero_boundary_flux_;
  size_t num_groups_;

public:
  explicit SweepBoundary(LBSBoundaryType bndry_type,
                         size_t num_groups,
                         CoordinateSystemType coord_type)
    : type_(bndry_type), coord_type_(coord_type), num_groups_(num_groups)
  {
    zero_boundary_flux_.resize(num_groups_, 0.0);
  }

  virtual ~SweepBoundary() = default;

  LBSBoundaryType Type() const { return type_; }

  CoordinateSystemType CoordType() const { return coord_type_; }

  bool IsReflecting() const { return type_ == LBSBoundaryType::REFLECTING; }

  double GetEvaluationTime() const { return evaluation_time_; }

  void SetEvaluationTime(double time) { evaluation_time_ = time; }

  /// Returns a pointer to the location of the incoming flux.
  virtual double* PsiIncoming(uint64_t cell_local_id,
                              unsigned int face_num,
                              unsigned int fi,
                              unsigned int angle_num,
                              int group_num,
                              size_t gs_ss_begin);

  /// Returns a pointer to the location of the outgoing flux.
  virtual double* PsiOutgoing(uint64_t cell_local_id,
                              unsigned int face_num,
                              unsigned int fi,
                              unsigned int angle_num,
                              size_t gs_ss_begin);

  virtual void UpdateAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss) {}

  virtual bool CheckAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss)
  {
    return true;
  }

  virtual void Setup(const MeshContinuum& grid, const AngularQuadrature& quadrature) {}

  double* ZeroFlux(int group_num) { return &zero_boundary_flux_[group_num]; }
};

/**
 * This boundary function class can be derived from to
 * provide a much more custom experience. This function
 * is called during Setup.
 */
class BoundaryFunction
{
public:
  /// Customized boundary function by calling a lua routine.
  virtual std::vector<double>
  Evaluate(size_t cell_global_id,
           int cell_material_id,
           unsigned int face_index,
           unsigned int face_node_index,
           const Vector3& face_node_location,
           const Vector3& face_node_normal,
           const std::vector<int>& quadrature_angle_indices,
           const std::vector<Vector3>& quadrature_angle_vectors,
           const std::vector<std::pair<double, double>>& quadrature_phi_theta_angles,
           const std::vector<int>& group_indices,
           double time) = 0;

  virtual ~BoundaryFunction() = default;
};

} // namespace opensn
