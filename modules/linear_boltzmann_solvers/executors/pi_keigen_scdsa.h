// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/executors/pi_keigen.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion.h"

namespace opensn
{
class VectorGhostCommunicator;

class PowerIterationKEigenSCDSA : public PowerIterationKEigen
{
private:
  int accel_pi_max_its_;
  double accel_pi_k_tol_;
  bool accel_pi_verbose_;
  std::shared_ptr<DiffusionSolver> diffusion_solver_ = nullptr;

  const std::string diffusion_solver_sdm_;
  std::shared_ptr<SpatialDiscretization> continuous_sdm_ptr_ = nullptr;
  bool requires_ghosts_ = false;
  struct GhostInfo
  {
    std::shared_ptr<VectorGhostCommunicator> vector_ghost_communicator = nullptr;
    std::map<int64_t, int64_t> ghost_global_id_2_local_map;
  };
  GhostInfo lbs_pwld_ghost_info_;

  double diff_accel_diffusion_l_abs_tol_;
  int diff_accel_diffusion_max_iters_;
  bool diff_accel_diffusion_verbose_;
  std::string diff_accel_diffusion_petsc_options_;

public:
  static InputParameters GetInputParameters();

  explicit PowerIterationKEigenSCDSA(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

  /// Copies only the scalar moments from an lbs primary flux moments vector.
  std::vector<double> CopyOnlyPhi0(const LBSGroupset& groupset, const std::vector<double>& phi_in);

  /// Copies back only the scalar moments to a lbs primary flux vector.
  void ProjectBackPhi0(const LBSGroupset& groupset,
                       const std::vector<double>& input,
                       std::vector<double>& output);

  /// Creates a ghost communicator and all associated information.
  GhostInfo MakePWLDVecGhostCommInfo(const SpatialDiscretization& sdm,
                                     const UnknownManager& uk_man);

  /**
   * This method takes an input vector that is the local version of a PWLD discrete space and then
   * makes it continuous by applying nodal averages.
   */
  static std::vector<double>
  NodallyAveragedPWLDVector(const std::vector<double>& input,
                            const SpatialDiscretization& pwld_sdm,
                            const SpatialDiscretization& pwlc_sdm,
                            const UnknownManager& uk_man,
                            const PowerIterationKEigenSCDSA::GhostInfo& ghost_info);
};

} // namespace opensn
