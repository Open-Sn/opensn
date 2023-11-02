#ifndef LBS_CURVILINEAR_SWEEPCHUNK_PWL2_H
#define LBS_CURVILINEAR_SWEEPCHUNK_PWL2_H

#include "B_DiscreteOrdinatesSolver/SweepChunks/AAH_SweepChunk.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

namespace lbs
{

/** A sweep-chunk in point-symmetric and axial-symmetric
 *  curvilinear coordinates. */
class SweepChunkPWLRZ : public lbs::AAH_SweepChunk
{
  //  Attributes
private:
  const std::vector<lbs::UnitCellMatrices>& secondary_unit_cell_matrices_;
  /** Unknown manager. */
  chi_math::UnknownManager unknown_manager_;
  /** Sweeping dependency angular intensity (for each polar level). */
  std::vector<double> psi_sweep_;
  /** Mapping from direction linear index to direction polar level. */
  std::map<unsigned int, unsigned int> map_polar_level_;
  /** Normal vector to determine symmetric boundary condition. */
  chi_mesh::Vector3 normal_vector_boundary_;

  // Runtime params
  const MatDbl*  Maux_ = nullptr;

  unsigned int polar_level_ = 0;
  double fac_diamond_difference_ = 0.0;
  double fac_streaming_operator_ = 0.0;

  //  Methods
public:
  /** Constructor. */
  SweepChunkPWLRZ(
    const chi_mesh::MeshContinuum& grid,
    const chi_math::SpatialDiscretization& discretization_primary,
    const std::vector<lbs::UnitCellMatrices>& unit_cell_matrices,
    const std::vector<lbs::UnitCellMatrices>& secondary_unit_cell_matrices,
    std::vector<lbs::CellLBSView>& cell_transport_views,
    std::vector<double>& destination_phi,
    std::vector<double>& destination_psi,
    const std::vector<double>& source_moments,
    lbs::LBSGroupset& groupset,
    const std::map<int, lbs::XSPtr>& xs,
    int num_moments,
    int max_num_cell_dofs);

protected:
  // operations
  void CellDataCallback();
  void DirectionDataCallback();
  void PostCellDirSweepCallback();

  // rz kernels
  void KernelFEMRZVolumetricGradientTerm();
  void KernelFEMRZUpwindSurfaceIntegrals();
};

} // namespace lbs

#endif // LBS_CURVILINEAR_SWEEPCHUNK_PWL2_H
