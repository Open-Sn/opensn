#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"

namespace opensn
{
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
  UnknownManager unknown_manager_;
  /** Sweeping dependency angular intensity (for each polar level). */
  std::vector<double> psi_sweep_;
  /** Mapping from direction linear index to direction polar level. */
  std::map<unsigned int, unsigned int> map_polar_level_;
  /** Normal vector to determine symmetric boundary condition. */
  Vector3 normal_vector_boundary_;

  // Runtime params
  const MatDbl* Maux_ = nullptr;

  unsigned int polar_level_ = 0;
  double fac_diamond_difference_ = 0.0;
  double fac_streaming_operator_ = 0.0;

  //  Methods
public:
  /** Constructor. */
  SweepChunkPWLRZ(const MeshContinuum& grid,
                  const SpatialDiscretization& discretization_primary,
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

  /**Cell data callback.*/
  void CellDataCallback();
  /**Direction data callback.*/
  void DirectionDataCallback();
  /**Applies diamond differencing on azimuthal directions.*/
  void PostCellDirSweepCallback();

  // rz kernels

  /**Assembles the volumetric gradient term.*/
  void KernelFEMRZVolumetricGradientTerm();
  /**Performs the integral over the surface of a face.*/
  void KernelFEMRZUpwindSurfaceIntegrals();
};

} // namespace lbs
} // namespace opensn
