// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/discrete_ordinates_keigen_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include <memory>
#include <limits>
#include <map>
#include <petscksp.h>
#include <string>
#include <tuple>
#include <vector>

namespace opensn
{

class CMFDAcceleration : public DiscreteOrdinatesKEigenAcceleration
{
public:
  explicit CMFDAcceleration(const InputParameters& params);
  ~CMFDAcceleration() override;

  void Initialize() final;
  void PreExecute() final;
  void PrePowerIteration() final;
  double PostPowerIteration() final;
  bool AllowsPowerIterationConvergence() const final { return last_update_allows_convergence_; }

  const CMFDCoarseMesh& GetCoarseMesh() const { return coarse_mesh_; }

  static InputParameters GetInputParameters();
  static std::shared_ptr<CMFDAcceleration> Create(const ParameterBlock& params);

private:
  struct BalanceResidual
  {
    double max_abs = 0.0;
    double relative_l2 = 0.0;
  };

  struct CoarseMeshDiagnostics
  {
    std::size_t local_fine_cells = 0;
    std::size_t local_coarse_cells = 0;
    std::size_t local_max_fine_cells_per_coarse_cell = 0;
    std::size_t local_max_faces_per_coarse_cell = 0;
    std::size_t local_undersized_coarse_cells = 0;
    double local_faces_per_coarse_cell_sum = 0.0;
    std::size_t global_fine_cells = 0;
    std::size_t global_coarse_cells = 0;
    std::size_t global_max_fine_cells_per_coarse_cell = 0;
    std::size_t global_max_faces_per_coarse_cell = 0;
    std::size_t global_undersized_coarse_cells = 0;
    double global_faces_per_coarse_cell_sum = 0.0;
  };

  struct FluxUpdateDiagnostics
  {
    double min_scalar_flux = std::numeric_limits<double>::max();
    bool has_nonfinite_flux = false;
    bool has_invalid_k = false;
  };

  void InitializeLinearSystem();
  void ApplyTransportUpdateScheme();
  void BuildCoarseFluxCache();
  void BuildFaceCurrentCache();
  void AssembleOperator();
  void AssembleRHS(Vec rhs, double k_eff, const std::vector<double>& coarse_source_phi) const;
  std::vector<double> SolveCoarseSystem(double k_eff,
                                        const std::vector<double>& coarse_source_phi) const;
  BalanceResidual ComputeBalanceResidual(double k_eff,
                                         const std::vector<double>& coarse_phi,
                                         const std::vector<double>& coarse_source_phi) const;
  double ComputeCoarseFissionProduction(const std::vector<double>& coarse_phi) const;
  std::size_t LocalUnknownCount() const;
  std::size_t GlobalUnknownCount() const;
  PetscInt MapDOF(uint64_t coarse_cell_global_id, unsigned int group) const;
  double ComputeOutwardCurrent(const CMFDCoarseCell& coarse_cell,
                               std::size_t face_index,
                               unsigned int group) const;
  void LogTimingSummary(double transport_time,
                        double restrict_old_time,
                        double restrict_new_time,
                        double coarse_flux_cache_time,
                        double face_current_cache_time,
                        double assemble_time,
                        double residual_time,
                        double fission_time,
                        double coarse_pi_time,
                        double update_time,
                        double post_time,
                        int coarse_pi_iterations,
                        int coarse_ksp_iterations,
                        const BalanceResidual& residual) const;
  CoarseMeshDiagnostics ComputeCoarseMeshDiagnostics() const;
  void LogCoarseMeshDiagnostics(const CoarseMeshDiagnostics& diagnostics) const;
  FluxUpdateDiagnostics AnalyzeFluxUpdate(const std::vector<double>& phi, double k_eff) const;
  bool IsAcceptableFluxUpdate(const FluxUpdateDiagnostics& diagnostics,
                              double min_allowed_scalar_flux) const;
  void UpdateAdaptiveRelaxation(double starting_relaxation,
                                double applied_damping,
                                bool skipped_correction);
  double ApplyFluxCorrectionWithDamping(const std::vector<double>& transport_phi_new,
                                        const std::vector<double>& coarse_phi,
                                        double transport_k_eff,
                                        double cmfd_k_eff,
                                        double& applied_damping,
                                        unsigned int& attempts,
                                        FluxUpdateDiagnostics& diagnostics,
                                        bool& skipped_correction);
  void LogCorrectionDiagnostics(double starting_relaxation,
                                double applied_damping,
                                unsigned int attempts,
                                const FluxUpdateDiagnostics& diagnostics,
                                bool skipped_correction) const;

  const std::string coarse_mesh_type_;
  const int aggregation_size_;
  const double relaxation_;
  const unsigned int correction_max_attempts_;
  const double correction_min_damping_;
  const double negative_flux_tolerance_;
  const bool adaptive_relaxation_;
  const double adaptive_relaxation_min_;
  const double adaptive_relaxation_max_;
  const double adaptive_relaxation_growth_;
  const double adaptive_relaxation_reduction_;
  const double adaptive_relaxation_accept_fraction_;
  const unsigned int adaptive_relaxation_successes_to_grow_;
  const unsigned int inactive_iterations_;
  const bool update_scheme_;
  const unsigned int update_wgs_max_its_;
  const double update_wgs_abs_tol_;
  const std::string coarse_solver_policy_;
  const std::size_t direct_coarse_solve_threshold_;
  const unsigned int first_group_ = 0;
  const unsigned int num_groups_;
  unsigned int outer_iteration_ = 0;
  bool last_update_allows_convergence_ = true;
  double current_relaxation_ = 0.0;
  unsigned int accepted_strong_corrections_ = 0;
  double last_restrict_old_time_ = 0.0;
  double transport_start_time_ = 0.0;
  CMFDCoarseMesh coarse_mesh_;
  std::vector<double> coarse_phi_old_;
  std::vector<double> coarse_phi_new_;
  std::map<std::tuple<uint64_t, unsigned int>, double> coarse_phi_cache_;
  std::map<std::tuple<uint64_t, uint64_t, unsigned int>, double> face_current_cache_;

  Mat A_ = nullptr;
  Vec rhs_ = nullptr;
  KSP ksp_ = nullptr;
  mutable int last_ksp_iterations_ = 0;
};

} // namespace opensn
