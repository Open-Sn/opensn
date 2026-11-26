// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/discrete_ordinates_keigen_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "framework/math/petsc_utils/petsc_utils.h"
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
class LogicalVolume;

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
  bool UsesSweepsPerUpdate() const final { return true; }

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
   * The matrix is condensed on demand from the cross-section objects of the last assembly, which
   * the problem shares. Replacing the problem's cross sections does not change it, but changing
   * those objects in place does.
   *
   * \param coarse_cell_global_id Global id of a coarse cell owned by the calling rank.
   * \return Row-major matrix whose entry [g, gp] is the transfer from coarse group gp to g. The
   *         diagonal is the within-group transfer, which the removal coefficient accounts for.
   * \throw std::logic_error If those objects no longer have the solver's group count.
   */
  std::vector<double> GetHomogenizedScatteringXS(uint64_t coarse_cell_global_id) const;
  /**
   * Return the production matrix of the last assembled CMFD operator for one coarse cell.
   *
   * The operator lifecycle and exceptions are as for GetHomogenizedRemovalXS(). Not collective.
   *
   * The matrix is condensed on demand, as for GetHomogenizedScatteringXS().
   *
   * \param coarse_cell_global_id Global id of a coarse cell owned by the calling rank.
   * \return Row-major matrix whose entry [g, gp] is the production in coarse group g per unit
   *         flux in coarse group gp, including delayed production when precursors are enabled.
   * \throw std::logic_error If those objects no longer have the solver's group count.
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

  /// Local counts of coarse operator entries whose sign is not that of an M-matrix.
  struct OperatorSignCounts
  {
    std::size_t nonpositive_diagonals = 0;
    std::size_t positive_off_diagonals = 0;
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
    unsigned int block_id = 0;
    double volume = 0.0;
    /// Per-transport-group weight of the region in the homogenization of its coarse group.
    std::vector<double> weights;
  };

  /// Homogenization data of one coarse cell. Transfer and production matrices are condensed on
  /// demand from the material weights so that storage stays O(materials x groups). The removal
  /// and diffusion vectors are indexed by coarse group.
  struct CoarseCellXS
  {
    std::vector<CoarseCellMaterial> materials;
    std::vector<double> removal;
    std::vector<double> diffusion;
  };

  /// The coarse mesh selected by the coarse_mesh parameter.
  struct CoarseMeshSpec
  {
    /// One of identity, local_aggregation, global_aggregation, reactor_lattice, cartesian, or
    /// logical_volumes.
    std::string type;
    /// Target fine cells per coarse cell of an aggregation.
    int size = 32;
    std::vector<double> x_edges;
    std::vector<double> y_edges;
    std::vector<double> z_edges;
    std::vector<std::shared_ptr<LogicalVolume>> logical_volumes;
    std::vector<double> lattice_origin;
    std::vector<double> lattice_pitch;
    std::vector<int> lattice_shape;
    std::vector<double> lattice_axial_edges;
  };

  /// Parses the coarse_mesh parameter and the deprecated aggregation_size. Throws on input that
  /// does not match the schema of the selected type.
  static CoarseMeshSpec ParseCoarseMeshSpec(const InputParameters& params);
  /// Number of fine groups per coarse group for the coarse_groups parameter and the deprecated
  /// group_aggregation_size.
  static unsigned int ParseGroupAggregationSize(const InputParameters& params,
                                                unsigned int num_groups);
  void InitializeLinearSystem();
  /// Chooses the coarse linear solver: LU while its estimated memory is small, FGMRES with an
  /// additive Schwarz preconditioner otherwise. Collective.
  void ConfigureCoarseLinearSolver();
  void ConfigureTransportSolve(unsigned int max_iterations, double residual_tolerance);
  /// Cache the diffusion coefficient and, once a transport iterate exists, the coarse flux of
  /// every local and face-neighbor coarse cell.
  void BuildCoarseFluxCache();
  /// The coarse flux a neighbor's closure sees: non-positive for a held value.
  double CachedCoarseFlux(std::size_t local_index) const;
  void BuildFaceCurrentCache();
  void InitializeHomogenization();
  /// Homogenize with the weights of phi, or with volume weights if use_flux is false. Collective.
  void BuildHomogenizedXS(const std::vector<double>& phi, bool use_flux);
  const CoarseCellXS& LocalHomogenizedXS(uint64_t coarse_cell_global_id) const;
  /// Throw if a cross-section object of the last assembly has changed its group count.
  void CheckAssembledGroupCount(const MultiGroupXS& xs, unsigned int block_id) const;
  void ComputeCoarseScattering(const CoarseCellXS& coarse_xs, std::vector<double>& matrix) const;
  /// Condense the production matrix; return false, with a zero matrix, if nothing is fissionable.
  bool ComputeCoarseProduction(const CoarseCellXS& coarse_xs, std::vector<double>& matrix) const;
  /**
   * Assemble the coarse operator without setting up the coarse linear solver. Before a transport
   * iterate exists the faces are closed by the diffusion coupling and the boundaries by the
   * Marshak condition. Collective.
   */
  OperatorSignCounts AssembleOperator(bool from_transport_iterate);
  /// Sets A_shift_ = A_ - F / k_shift and makes it the operator of the coarse solves. Collective.
  void SetShiftedOperator(double k_shift);
  /**
   * Warns once if coarse cells are optically thick. CMFD with one sweep per update converges at
   * a rate that approaches the scattering ratio as the coarse cells thicken. Collective.
   */
  void CheckOpticalThickness();
  /// Production right-hand side, plus the held values and their sources times `held_scale`.
  void AssembleRHS(Vec rhs,
                   double k_eff,
                   const std::vector<double>& coarse_source_phi,
                   double held_scale) const;
  std::vector<double> SolveCoarseSystem(double k_eff,
                                        const std::vector<double>& coarse_source_phi,
                                        double held_scale) const;
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
  /**
   * Prolongs the coarse solution onto the transport iterate. A coarse value whose ratio would
   * change the sign of the fine flux keeps the transport iterate, scaled by `coarse_amplitude`,
   * the ratio of the coarse solution's amplitude to the transport iterate's. A correction that
   * gives a non-finite flux or an invalid k is skipped and the transport iterate is kept. Returns
   * the k of the iterate that is kept. Collective.
   */
  double ApplyFluxCorrection(const std::vector<double>& transport_phi_new,
                             const std::vector<double>& coarse_phi,
                             double coarse_amplitude,
                             double old_fission_production,
                             double transport_k_eff,
                             FluxUpdateDiagnostics& diagnostics,
                             bool& skipped_correction,
                             std::string& skip_reason);

  const CoarseMeshSpec coarse_mesh_spec_;
  /// The deprecated update_wgs_max_its when given, which then overrides the sweeps of the owning
  /// solver's automatic choice.
  const unsigned int legacy_update_wgs_max_its_;
  const unsigned int first_group_ = 0;
  const unsigned int num_groups_;
  const unsigned int group_aggregation_size_;
  const unsigned int num_coarse_groups_;
  unsigned int outer_iteration_ = 0;
  /// Transport sweeps per update, from the owning solver at each Execute().
  unsigned int sweeps_per_update_ = 1;
  bool last_update_allows_convergence_ = true;
  bool last_balance_residual_known_ = false;
  bool last_correction_skipped_ = false;
  unsigned int num_applied_corrections_ = 0;
  unsigned int num_skipped_corrections_ = 0;
  unsigned int consecutive_skipped_corrections_ = 0;
  bool warned_persistent_skips_ = false;
  bool warned_homogenization_fallback_ = false;
  bool warned_coarse_pi_limit_ = false;
  bool warned_kept_coarse_values_ = false;
  /// Per local coarse value, 1 if the coarse solve holds it at its restricted transport flux (see
  /// AssembleOperator).
  std::vector<char> held_;
  /// Per local coarse value: a held value's restricted transport flux, or a free value's incoming
  /// partial currents from held neighbors.
  std::vector<double> held_source_;
  /// Global number of coarse values that kept the transport flux in the last correction.
  std::size_t last_kept_coarse_values_ = 0;
  bool checked_optical_thickness_ = false;
  /// Global count of face closures that fell back to the diffusion coupling alone because a
  /// coarse flux was zero, over the last solve.
  std::size_t num_closure_fallbacks_ = 0;
  std::string last_correction_skip_reason_;
  double last_transport_current_balance_relative_l2_ = std::numeric_limits<double>::infinity();
  double last_restrict_old_time_ = 0.0;
  double transport_start_time_ = 0.0;
  double old_fission_production_ = 0.0;
  bool using_direct_coarse_solver_ = false;
  /// Whether petsc_options configure the coarse linear solver instead of the automatic choice.
  bool coarse_solver_from_options_ = false;
  CMFDCoarseMesh coarse_mesh_;
  std::vector<double> coarse_phi_old_;
  std::vector<double> coarse_phi_new_;
  std::map<std::tuple<uint64_t, unsigned int>, double> coarse_phi_cache_;
  /// Homogenization data of the last assembled operator, per local coarse cell.
  std::vector<CoarseCellXS> homogenized_xs_;
  /// Cross sections and mode of the last assembled operator.
  BlockID2XSMap assembled_xs_map_;
  bool assembled_adjoint_ = false;
  bool assembled_delayed_production_ = false;
  /// Per local fine cell, the destination rank and record of its homogenization flux integrals.
  std::vector<std::pair<int, std::size_t>> homogenization_send_records_;
  /// Number of (coarse cell, block) homogenization records sent to each rank.
  std::map<int, std::size_t> homogenization_send_counts_;
  /// Per source rank, the (coarse local id, material index) of each received record.
  std::map<int, std::vector<std::pair<std::size_t, std::size_t>>> homogenization_receive_targets_;
  std::map<std::tuple<uint64_t, unsigned int>, double> diffusion_coefficient_cache_;
  std::map<std::tuple<uint64_t, std::size_t, uint64_t, unsigned int>, double> face_current_cache_;

  Mat A_ = nullptr;
  /// A_ minus the coarse production matrix over the Wielandt shift, for the coarse solves.
  Mat A_shift_ = nullptr;
  Vec rhs_ = nullptr;
  KSP ksp_ = nullptr;
  /// Options of the coarse solver, which are in the global database only while it uses them.
  std::unique_ptr<PETScOptionsDatabase> petsc_options_db_;
  mutable int last_ksp_iterations_ = 0;
  mutable bool last_ksp_converged_ = true;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<CMFDAcceleration> Create(const ParameterBlock& params);
};

} // namespace opensn
