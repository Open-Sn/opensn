#pragma once

#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/sweep_boundary/sweep_boundary.h"
#include <vector>
#include <limits>

namespace opensn
{

/**
 * Specified incident fluxes homogenous on a boundary.
 */
class BoundaryIncidentHeterogeneous : public SweepBoundary
{
private:
  std::unique_ptr<BoundaryFunction> boundary_function_;
  const uint64_t ref_boundary_id_;

  typedef std::vector<double> FaceNodeData;
  typedef std::vector<FaceNodeData> FaceData;
  typedef std::vector<FaceData> CellData;

  std::vector<CellData> local_cell_data_;

public:
  explicit BoundaryIncidentHeterogeneous(
    size_t in_num_groups,
    std::unique_ptr<BoundaryFunction> in_bndry_function,
    uint64_t in_ref_boundary_id,
    CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::INCIDENT_ANISOTROPIC_HETEROGENEOUS, in_num_groups, coord_type),
      boundary_function_(std::move(in_bndry_function)),
      ref_boundary_id_(in_ref_boundary_id)
  {
  }

  double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   int group_num,
                                   size_t gs_ss_begin) override;

  void Setup(const MeshContinuum& grid, const AngularQuadrature& quadrature) override;
};

} // namespace opensn
