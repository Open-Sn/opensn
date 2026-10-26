// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/discrete_ordinates_keigen_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include <petscksp.h>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace opensn
{

class MultiGroupXS;

class CMFDAcceleration : public DiscreteOrdinatesKEigenAcceleration
{
public:
  explicit CMFDAcceleration(const InputParameters& params);
  ~CMFDAcceleration() override;

  void CheckRequirements(std::vector<std::string>& errors) const override;
  void Initialize() final;
  void PreExecute() final;
  void PrePowerIteration() final;
  double PostPowerIteration() final;
  bool AllowsPowerIterationConvergence() const final { return last_update_allows_convergence_; }
  std::string GetPowerIterationConvergenceInfo() const final;

  const CMFDCoarseMesh& GetCoarseMesh() const { return coarse_mesh_; }
  /// Global number of CMFD coarse cells; zero before Initialize().
  std::size_t GetNumCoarseCells() const { return coarse_mesh_.NumGlobalCells(); }
  /**
   * Return the removal coefficients of the last assembled CMFD operator for one coarse cell.
   *
   * Initialize(), which the owning solver calls at the start of each Execute(), assembles a
   * volume-weighted operator; each active power iteration reassembles it from that iteration's
   * transport iterate before the correction is applied. Not collective.
   *
   * \param coarse_cell_global_id Global id of a coarse cell owned by the calling rank.
   * \return One removal coefficient per CMFD energy group.
   * \throw std::invalid_argument If the coarse cell is not owned by the calling rank.
   * \throw std::logic_error If called before Initialize().
   */
  std::vector<double> GetHomogenizedRemovalXS(uint64_t coarse_cell_global_id) const;
  /**
   * Return the diffusion coefficients of the last assembled CMFD operator for one coarse cell.
   *
   * The operator lifecycle and exceptions are as for GetHomogenizedRemovalXS(). Not collective.
   *
   * \param coarse_cell_global_id Global id of a coarse cell owned by the calling rank.
   * \return One diffusion coefficient per CMFD energy group.
   */
  std::vector<double> GetHomogenizedDiffusionCoefficient(uint64_t coarse_cell_global_id) const;
  /**
   * Return the P0 transfer matrix of the last assembled CMFD operator for one coarse cell.
   *
   * The operator lifecycle and exceptions are as for GetHomogenizedRemovalXS(). Not collective.
   *
   * \param coarse_cell_global_id Global id of a coarse cell owned by the calling rank.
   * \return Row-major matrix whose entry [g, gp] is the transfer from coarse group gp to g. The
   *         diagonal is the within-group transfer, which the removal coefficient accounts for.
   */
  std::vector<double> GetHomogenizedScatteringXS(uint64_t coarse_cell_global_id) const;
  /**
   * Return the production matrix of the last assembled CMFD operator for one coarse cell.
   *
   * The operator lifecycle and exceptions are as for GetHomogenizedRemovalXS(). Not collective.
   *
   * \param coarse_cell_global_id Global id of a coarse cell owned by the calling rank.
   * \return Row-major matrix whose entry [g, gp] is the production in coarse group g per unit
   *         flux in coarse group gp, including delayed production when precursors are enabled.
   */
  std::vector<double> GetHomogenizedProductionXS(uint64_t coarse_cell_global_id) const;
  /// Number of active CMFD corrections applied during the last power-iteration solve.
  unsigned int GetNumAppliedCorrections() const { return num_applied_corrections_; }
  /// Number of active CMFD corrections skipped during the last power-iteration solve.
  unsigned int GetNumSkippedCorrections() const { return num_skipped_corrections_; }

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

  /// The fine cells of one material block within a coarse cell.
  struct CoarseCellMaterial
  {
    /// Material block id of the region.
    unsigned int block_id = 0;
    /// Physical volume of the region.
    double volume = 0.0;
    /// Per-transport-group weight of the region in the homogenization of its coarse group.
    std::vector<double> weights;
  };

  /// Homogenization data of one coarse cell. Transfer and production matrices are condensed on
  /// demand from the material weights so that storage stays O(materials x groups).
  struct CoarseCellXS
  {
    /// Material regions of the coarse cell.
    std::vector<CoarseCellMaterial> materials;
    /// Removal coefficient per coarse group.
    std::vector<double> removal;
    /// Diffusion coefficient per coarse group.
    std::vector<double> diffusion;
  };

  bool EvaluateBalanceConvergence(const BalanceResidual& transport_current_balance_residual,
                                  bool skipped_correction) const;
  void InitializeLinearSystem();
  void ConfigureCoarseLinearSolver(bool direct);
  void ConfigureTransportSolve(unsigned int max_iterations, double residual_tolerance);
  /// Cache the diffusion coefficient and, once a transport iterate exists, the coarse flux of
  /// every local and face-neighbor coarse cell.
  void BuildCoarseFluxCache();
  void BuildFaceCurrentCache();
  void InitializeHomogenization();
  void BuildHomogenizedXS(const std::vector<double>& phi);
  const CoarseCellXS& LocalHomogenizedXS(uint64_t coarse_cell_global_id) const;
  void ComputeCoarseScattering(const CoarseCellXS& coarse_xs, std::vector<double>& matrix) const;
  /// Condense the production matrix; return false, with a zero matrix, if nothing is fissionable.
  bool ComputeCoarseProduction(const CoarseCellXS& coarse_xs, std::vector<double>& matrix) const;
  void AssembleOperator();
  void AssembleRHS(Vec rhs, double k_eff, const std::vector<double>& coarse_source_phi) const;
  std::vector<double> SolveCoarseSystem(double k_eff,
                                        const std::vector<double>& coarse_source_phi) const;
  BalanceResidual ComputeBalanceResidual(double k_eff,
                                         const std::vector<double>& coarse_phi,
                                         const std::vector<double>& coarse_source_phi) const;
  BalanceResidual
  ComputeTransportCurrentBalanceResidual(double k_eff,
                                         const std::vector<double>& coarse_phi,
                                         const std::vector<double>& coarse_source_phi) const;
  double ComputeCoarseFissionProduction(const std::vector<double>& coarse_phi) const;
  std::size_t LocalUnknownCount() const;
  std::size_t GlobalUnknownCount() const;
  PetscInt MapDOF(uint64_t coarse_cell_global_id, unsigned int coarse_group) const;
  double ComputeOutwardCurrent(const CMFDCoarseCell& coarse_cell,
                               std::size_t face_index,
                               unsigned int coarse_group) const;
  std::pair<double, double> ComputePartialOutwardCurrents(const CMFDCoarseCell& coarse_cell,
                                                          std::size_t face_index,
                                                          unsigned int coarse_group) const;
  unsigned int FineGroupBegin(unsigned int coarse_group) const;
  unsigned int FineGroupEnd(unsigned int coarse_group) const;
  FluxUpdateDiagnostics AnalyzeFluxUpdate(const std::vector<double>& phi, double k_eff) const;
  double ApplyFluxCorrectionWithDamping(const std::vector<double>& transport_phi_new,
                                        const std::vector<double>& coarse_phi,
                                        double old_fission_production,
                                        double transport_k_eff,
                                        double& applied_damping,
                                        unsigned int& attempts,
                                        FluxUpdateDiagnostics& diagnostics,
                                        bool& skipped_correction,
                                        std::string& skip_reason);

  const std::string coarse_mesh_type_;
  const std::string current_closure_;
  const int aggregation_size_;
  /// Whether aggregated coarse cells may cross material-block boundaries.
  const bool allow_mixed_material_aggregation_;
  const unsigned int group_aggregation_size_;
  const double relaxation_;
  const unsigned int correction_max_attempts_;
  const double correction_min_damping_;
  const double negative_flux_tolerance_;
  const unsigned int inactive_iterations_;
  const unsigned int update_wgs_max_its_;
  const double update_wgs_abs_tol_;
  const double balance_residual_tolerance_;
  const double l_rel_tol_;
  const std::string coarse_solver_policy_;
  const std::size_t coarse_direct_solve_threshold_;
  const unsigned int first_group_ = 0;
  const unsigned int num_groups_;
  const unsigned int num_coarse_groups_;
  unsigned int outer_iteration_ = 0;
  bool last_update_allows_convergence_ = true;
  bool last_balance_residual_known_ = false;
  bool last_correction_skipped_ = false;
  unsigned int num_applied_corrections_ = 0;
  unsigned int num_skipped_corrections_ = 0;
  unsigned int consecutive_skipped_corrections_ = 0;
  bool warned_persistent_skips_ = false;
  /// Whether the volume-weight homogenization fallback was reported in this Execute().
  bool warned_homogenization_fallback_ = false;
  std::string last_correction_skip_reason_;
  double last_transport_current_balance_relative_l2_ = std::numeric_limits<double>::infinity();
  double last_restrict_old_time_ = 0.0;
  double transport_start_time_ = 0.0;
  double old_fission_production_ = 0.0;
  bool using_direct_coarse_solver_ = false;
  CMFDCoarseMesh coarse_mesh_;
  std::vector<double> coarse_phi_old_;
  std::vector<double> coarse_phi_new_;
  std::map<std::tuple<uint64_t, unsigned int>, double> coarse_phi_cache_;
  /// Homogenization data of the last assembled operator, per local coarse cell.
  std::vector<CoarseCellXS> homogenized_xs_;
  /// Per local fine cell, the destination rank and record of its homogenization flux integrals.
  std::vector<std::pair<int, std::size_t>> homogenization_send_records_;
  /// Number of (coarse cell, block) homogenization records sent to each rank.
  std::map<int, std::size_t> homogenization_send_counts_;
  /// Per source rank, the (coarse local id, material index) of each received record.
  std::map<int, std::vector<std::pair<std::size_t, std::size_t>>> homogenization_receive_targets_;
  /// Diffusion coefficient per (coarse cell, coarse group) for local and neighbor coarse cells.
  std::map<std::tuple<uint64_t, unsigned int>, double> diffusion_coefficient_cache_;
  std::map<std::tuple<uint64_t, std::size_t, uint64_t, unsigned int>, double> face_current_cache_;

  Mat A_ = nullptr;
  Vec rhs_ = nullptr;
  KSP ksp_ = nullptr;
  mutable int last_ksp_iterations_ = 0;
  mutable bool last_ksp_converged_ = true;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<CMFDAcceleration> Create(const ParameterBlock& params);
};

} // namespace opensn
