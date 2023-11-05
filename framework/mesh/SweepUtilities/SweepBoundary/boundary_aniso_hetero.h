#pragma once

#include "mesh/chi_mesh.h"
#include "math/chi_math.h"
#include "sweep_boundary.h"
#include <vector>
#include <limits>

namespace chi_mesh::sweep_management
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
    chi_math::CoordinateSystemType coord_type = chi_math::CoordinateSystemType::CARTESIAN)
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

  void Setup(const chi_mesh::MeshContinuum& grid,
             const chi_math::AngularQuadrature& quadrature) override;
};

} // namespace chi_mesh::sweep_management
