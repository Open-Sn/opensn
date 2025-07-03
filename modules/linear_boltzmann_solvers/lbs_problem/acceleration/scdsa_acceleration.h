// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/lbs_keigen_acceleration.h"

namespace opensn
{

class DiffusionSolver;
class SpatialDiscretization;
class LinearSolver;
struct WGSContext;

class SCDSAAcceleration : public LBSKEigenAcceleration
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
  /// Spatial discretization method, from parameters
  const std::string sdm_;

  /// Front groupset WGS solver from the LBSProblem
  std::shared_ptr<LinearSolver> front_wgs_solver_;
  /// Front groupset WGSContext from the LBSProblem
  std::shared_ptr<WGSContext> front_wgs_context_;

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
  ///@}
};
} // namespace opensn
