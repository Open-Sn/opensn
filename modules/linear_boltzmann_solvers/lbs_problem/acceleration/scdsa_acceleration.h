// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/lbs_acceleration.h"

namespace opensn
{

class DiffusionSolver;
class SpatialDiscretization;
class VectorGhostCommunicator;
class LinearSolver;
struct WGSContext;

class SCDSAAcceleration : public LBSAcceleration
{
public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<SCDSAAcceleration> Create(const ParameterBlock& params);

  explicit SCDSAAcceleration(const InputParameters& params);

  void Initialize() override final;
  void PreExecute() override final;
  void PrePowerIteration() override final;
  double PostPowerIteration() override final;

private:
  struct GhostInfo
  {
    std::shared_ptr<VectorGhostCommunicator> vector_ghost_communicator = nullptr;
    std::map<int64_t, int64_t> ghost_global_id_2_local_map;
  };

  /// Creates a ghost communicator and all associated information.
  GhostInfo MakePWLDVecGhostCommInfo() const;

  /// Copies only the scalar moments from an lbs primary flux moments vector.
  void CopyOnlyPhi0(const std::vector<double>& phi_in, std::vector<double>& phi_local);

  /// Copies back only the scalar moments to a lbs primary flux vector.
  void ProjectBackPhi0(const std::vector<double>& input, std::vector<double>& output) const;

  /**
   * This method takes an input vector that is the local version of a PWLD discrete space and then
   * makes it continuous by applying nodal averages.
   */
  void NodallyAveragedPWLDVector(const std::vector<double>& input,
                                 std::vector<double>& output) const;

  /// Spatial discretization method, from parameters
  const std::string sdm_;
  /// Maximum inner power iteration count, from parameters
  const int pi_max_its_;
  /// k-eigenvalue tolerance for inner power iterations, from parameters
  const double pi_k_tol_;

  /// Front groupset from the LBSProblem
  LBSGroupset& front_gs_;
  /// Front groupset WGS solver from the LBSProblem
  std::shared_ptr<LinearSolver> front_wgs_solver_;
  /// Front groupset WGSContext from the LBSProblem
  std::shared_ptr<WGSContext> front_wgs_context_;

  /// Underlying diffusion solver, set during Initialize()
  std::shared_ptr<DiffusionSolver> diffusion_solver_ = nullptr;
  /// Underlying continuous discretization (for pwlc), set during Initialize()
  std::shared_ptr<SpatialDiscretization> continuous_sdm_ptr_ = nullptr;
  /// Ghost information, needed with a continuous discretization, set during Initialize()
  GhostInfo lbs_pwld_ghost_info_;

  /// Whether or not an extra sweep is required, set during PreExecute()
  bool extra_sweep_ = false;

  /**
   * Work vectors
   */
  ///@{
  std::vector<double> phi_temp_;
  std::vector<double> Sf_ell_;
  std::vector<double> Sf0_ell_;
  std::vector<double> phi0_ell_;
  std::vector<double> phi0_star_;
  std::vector<double> Ss0_res_;
  std::vector<double> epsilon_k_;
  std::vector<double> epsilon_kp1_;
  std::vector<double> Sf0_aux_;
  std::vector<double> Ss0_;
  std::vector<double> copy_only_phi0_tmp_;
  ///@}
};
} // namespace opensn
