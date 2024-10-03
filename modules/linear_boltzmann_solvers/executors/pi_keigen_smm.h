// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/ndarray.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion.h"

namespace opensn
{
class VectorGhostCommunicator;
class GhostedParallelSTLVector;

class PowerIterationKEigenSMM : public PowerIterationKEigen
{
protected:
  struct GhostInfo
  {
    std::shared_ptr<VectorGhostCommunicator> vector_ghost_communicator;
    std::map<int64_t, int64_t> ghost_global_to_local_map;
  };

public:
  static InputParameters GetInputParameters();
  explicit PowerIterationKEigenSMM(const InputParameters& params);

  void Initialize() override;
  void Execute() override;

protected:
  void ComputeClosures(const std::vector<std::vector<double>>& psi);
  std::vector<double> ComputeSourceCorrection() const;

  void AssembleDiffusionBCs() const;
  std::vector<double> AssembleDiffusionRHS(const std::vector<double>& q0) const;
  std::vector<double> SetNodalDiffusionFissionSource(const std::vector<double>& phi0) const;
  std::vector<double> SetNodalDiffusionScatterSource(const std::vector<double>& phi0) const;

  void ComputeAuxiliaryUnitCellMatrices();
  void ComputeBoundaryFactors();

  /**
   * Transfer a transport flux moments vector to a diffusion vector.
   *
   * \param[in] input A transport flux moments vector to transfer data from.
   * \param[out] output A vector to place the transferred data.
   */
  void TransferTransportToDiffusion(const std::vector<double>& input,
                                    std::vector<double>& output) const;

  /**
   * Transfer a diffusion scalar flux vector to a transport flux moments vector.
   *
   * \param[in] input A diffusion scalar flux vector to transfer data from.
   * \param[out] output A transport flux moments vector to transfer data to.
   */
  void TransferDiffusionToTransport(const std::vector<double>& input,
                                    std::vector<double>& output) const;

  double CheckScalarFluxConvergence(const std::vector<double>& phi_new,
                                    const std::vector<double>& phi_old);

protected:
  unsigned int dimension_;
  std::vector<std::vector<double>>& psi_new_local_;

  // Second moment closures
  UnknownManager tensor_uk_man_;
  std::shared_ptr<GhostedParallelSTLVector> tensors_;
  std::map<uint64_t, std::vector<double>> betas_;

  // Quadrature approximated boundary factors per groupset
  std::map<uint64_t, std::vector<double>> bndry_factors_;

  // Cell-wise tensor stiffness matrices.
  std::vector<NDArray<double, 4>> K_tensor_matrices_;

  // Diffusion solver attributes
  std::shared_ptr<DiffusionSolver> diffusion_solver_ = nullptr;
  std::shared_ptr<SpatialDiscretization> pwlc_ptr_ = nullptr;

  bool ghosts_required_;
  GhostInfo ghost_info_;

  // Additional options
  const unsigned int accel_pi_max_its_;
  const double accel_pi_k_tol_;
  const bool accel_pi_verbose_;

  const std::string diffusion_sdm_name_;
  const int diffusion_l_max_its_;
  const double diffusion_l_abs_tol_;
  const std::string diffusion_petsc_options_;
  const bool diffusion_verbose_;

protected:
  /**
   * Return a PWLD vector whose entries contain the average value of
   * all associated discontinuous nodes.
   *
   * \param pwld_vector A vector on a PWLD discretization.
   * \param pwld The PWLD discretization.
   * \param pwlc The PWLC discretization.
   * \param uk_man The unknown manager for the PWLD unknowns.
   * \param ghost_info The ghost node info for the PWLC unknowns.
   */
  static std::vector<double>
  ComputeNodallyAveragedPWLDVector(const std::vector<double>& pwld_vector,
                                   const SpatialDiscretization& pwld,
                                   const SpatialDiscretization& pwlc,
                                   const UnknownManager& uk_man,
                                   const GhostInfo& ghost_info);

  static std::vector<int64_t> MakePWLDGhostIndices(const SpatialDiscretization& pwld,
                                                   const UnknownManager& uk_man);

  /**
   * Create a vector ghost communicator and a global to local mapping for
   * PWLD ghost cell data.
   *
   * \param pwld The PWLD discretization to develop the mapping for.
   * \param uk_man An unknown manager to develop the mapping for.
   */
  static GhostInfo MakePWLDGhostInfo(const SpatialDiscretization& pwld,
                                     const UnknownManager& uk_man);

  /**
   * Obtain the index of the associated node on the opposing side of a face.
   * Given a node, go through the nodes of a neighbor cell to find a matching
   * node. This routine returns the local cell node index of the specified node
   * for the neighboring cell.
   *
   * \param node The node to find an associated node for.
   * \param nbr_nodes The cell nodes on a neighbor cell.
   * \param epsilon A matching tolerance.
   * \return The local cell node index on the neighboring cell.
   */
  static int MapAssociatedFaceNode(const Vector3& node,
                                   const std::vector<Vector3>& nbr_nodes,
                                   double epsilon = 1.0e-12);
};

} // namespace opensn
