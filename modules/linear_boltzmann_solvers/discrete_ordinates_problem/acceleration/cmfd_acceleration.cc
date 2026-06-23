// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_vector_tools.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/linear_solver/petsc_linear_system_solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"
#include "framework/utils/timer.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/classic_richardson.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <set>
#include <stdexcept>
#include <utility>

namespace opensn
{
namespace
{

double
ProjectedDistanceToFace(const Vector3& cell_centroid,
                        const Vector3& face_centroid,
                        const Vector3& face_normal)
{
  return std::max(1.0e-12, std::fabs((face_centroid - cell_centroid).Dot(face_normal)));
}

double
InterfaceDiffusionCoefficient(const double owner_diffusion_coefficient,
                              const double neighbor_diffusion_coefficient,
                              const double owner_distance_to_face,
                              const double neighbor_distance_to_face,
                              const double face_area)
{
  const double denominator = owner_distance_to_face / owner_diffusion_coefficient +
                             neighbor_distance_to_face / neighbor_diffusion_coefficient;
  return face_area / std::max(1.0e-12, denominator);
}

double
NowSeconds()
{
  return program_timer.GetTime() * 1.0e-3;
}

double
MaxAcrossRanks(const double local_value)
{
  double global_value = 0.0;
  opensn::mpi_comm.all_reduce(local_value, global_value, mpi::op::max<double>());
  return global_value;
}

void
LogCMFDAcceleratorMetric(const std::string& category,
                         const int iteration,
                         const std::string& metric,
                         const std::string& value)
{
  std::ostringstream out;
  out << program_timer.GetTimeString() << " CMFD_ACCEL category=" << category;
  if (iteration >= 0)
    out << " iteration=" << iteration;
  out << " metric=" << metric << " value=" << value;
  log.Log() << no_wrap << out.str();
}

std::string
MetricValue(const double value)
{
  std::ostringstream out;
  out << std::scientific << std::setprecision(6) << value;
  return out.str();
}

template <typename T>
std::string
MetricValue(const T value)
{
  std::ostringstream out;
  out << value;
  return out.str();
}

double
RelativeL2Difference(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
  OpenSnInvalidArgumentIf(lhs.size() != rhs.size(), "RelativeL2Difference vector size mismatch.");
  double numerator = 0.0;
  double denominator = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i)
  {
    const double difference = lhs[i] - rhs[i];
    numerator += difference * difference;
    denominator += rhs[i] * rhs[i];
  }
  return std::sqrt(numerator / std::max(denominator, 1.0e-300));
}

} // namespace

OpenSnRegisterObjectInNamespace(lbs, CMFDAcceleration);

InputParameters
CMFDAcceleration::GetInputParameters()
{
  auto params = DiscreteOrdinatesKEigenAcceleration::GetInputParameters();

  params.AddOptionalParameter(
    "coarse_mesh",
    "local_aggregation",
    "CMFD coarse-mesh construction method. local_aggregation builds rank-local aggregates. "
    "global_aggregation builds connected same-block aggregates that may span MPI ranks. identity "
    "creates one CMFD cell per transport cell and is intended mainly for debugging and method "
    "comparisons.");
  params.AddOptionalParameter("current_closure",
                              "auto",
                              "CMFD face-current closure. Valid choices are auto, net, and "
                              "partial. auto probes the early CMFD operator behavior and may use "
                              "net, partial, or a blend; fixed choices are useful for comparison "
                              "studies or cases where auto is not robust.");
  params.AddOptionalParameter(
    "aggregation_size",
    32,
    "Target number of fine cells per rank-local aggregated CMFD coarse cell. Larger values "
    "reduce CMFD cost but coarsen the spatial correction; smaller values increase coarse-solve "
    "cost but can improve robustness.");
  params.AddOptionalParameter(
    "group_aggregation_size",
    1,
    "Number of transport energy groups per CMFD coarse energy group, not the final number of "
    "coarse groups. Use (num_groups + target_coarse_groups - 1) / target_coarse_groups to target "
    "a desired number of CMFD energy groups.");
  params.AddOptionalParameter(
    "relaxation",
    0.5,
    "Requested relaxation factor for the CMFD scalar-flux correction. The correction limiter may "
    "damp or skip an individual correction if the requested correction is not admissible.");
  params.AddOptionalParameter(
    "correction_max_attempts",
    10,
    "[developer/debug] Maximum CMFD correction damping attempts before skipping the "
    "correction for the current transport update. Each failed attempt halves the damping.");
  params.AddOptionalParameter(
    "correction_min_damping",
    1.0e-4,
    "[developer/debug] Minimum CMFD correction damping factor considered during correction "
    "limiting. If no admissible correction is found above this damping, the CMFD correction is "
    "skipped for the current transport update.");
  params.AddOptionalParameter("negative_flux_tolerance",
                              1.0e-6,
                              "[developer/debug] Allowed scalar-flux undershoot for accepting a "
                              "CMFD correction. Corrections that drive scalar flux below this "
                              "guard are damped or skipped.");
  params.AddOptionalParameter(
    "inactive_iterations",
    0,
    "[developer/debug] Number of initial power iterations before applying CMFD "
    "corrections. Transport update controls are still active during "
    "inactive iterations.");
  params.AddOptionalParameter(
    "update_wgs_max_its",
    1,
    "Maximum WGS iterations used for each transport update before applying the CMFD "
    "correction. The default is one transport update iteration per power iteration.");
  params.AddOptionalParameter(
    "update_wgs_abs_tol",
    1.0e-12,
    "WGS absolute tolerance used for each transport update. When update_wgs_max_its=1, "
    "this tolerance is normally not the stopping criterion because only one WGS iteration "
    "is allowed.");
  params.AddOptionalParameter(
    "balance_residual_tolerance",
    1.0e-6,
    "Restricted transport-current balance residual tolerance required before "
    "CMFD allows power-iteration convergence. This is a second convergence "
    "guard in addition to the outer k_eff tolerance, not an estimate of the "
    "k-eigenvalue error. Large 3D problems may need a looser value than the "
    "outer k_eff tolerance.");
  params.AddOptionalParameter(
    "coarse_solver_policy",
    "auto",
    "[developer/debug] CMFD coarse solver policy. auto chooses direct PETSc LU below "
    "coarse_direct_solve_threshold global unknowns and iterative GMRES/Jacobi otherwise. "
    "direct, iterative, and petsc_options force a specific path. CMFD corrections from "
    "unconverged coarse linear solves are always skipped.");
  params.AddOptionalParameter("coarse_direct_solve_threshold",
                              20000,
                              "[developer/debug] Maximum global CMFD unknown count for using the "
                              "direct PETSc LU path when coarse_solver_policy is auto. Larger "
                              "values make auto use direct solves for larger coarse systems; "
                              "choose values based on available host memory and acceptable "
                              "factorization cost.");
  params.ConstrainParameterRange(
    "coarse_mesh",
    AllowableRangeList::New({"identity", "local_aggregation", "global_aggregation"}));
  params.ConstrainParameterRange("current_closure",
                                 AllowableRangeList::New({"auto", "net", "partial"}));
  params.ConstrainParameterRange(
    "coarse_solver_policy",
    AllowableRangeList::New({"auto", "direct", "iterative", "petsc_options"}));
  params.ConstrainParameterRange("aggregation_size", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("group_aggregation_size", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("relaxation", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("correction_max_attempts", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("correction_min_damping", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("negative_flux_tolerance", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("inactive_iterations", AllowableRangeLowLimit::New(0));
  params.ConstrainParameterRange("update_wgs_max_its", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("update_wgs_abs_tol", AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("balance_residual_tolerance", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("coarse_direct_solve_threshold", AllowableRangeLowLimit::New(1));
  params.ChangeExistingParamToOptional("name", "CMFDAcceleration");
  params.ChangeExistingParamToOptional(
    "l_abs_tol", 1.0e-7, "[developer/debug] Absolute residual tolerance for CMFD coarse solves.");
  params.ChangeExistingParamToOptional(
    "max_iters",
    100,
    "[developer/debug] Maximum iterations for each CMFD coarse linear solve. If a coarse "
    "linear solve does not converge within this limit, CMFD skips the correction rather than "
    "applying an untrusted low-order update. If this happens with an iterative coarse solve, "
    "use direct for modest systems, increase this limit, or provide stronger PETSc options.");
  params.ChangeExistingParamToOptional(
    "verbose", false, "[developer/debug] If true, prints detailed CMFD diagnostics and timings.");
  params.ChangeExistingParamToOptional(
    "petsc_options", std::string(""), "[developer/debug] Additional PETSc options.");
  params.ChangeExistingParamToOptional(
    "pi_max_its",
    50,
    "[developer/debug] Maximum inner power iterations for the CMFD coarse k solve.");
  params.ChangeExistingParamToOptional(
    "pi_k_tol",
    1.0e-8,
    "[developer/debug] k-eigenvalue tolerance for CMFD coarse power iterations.");

  return params;
}

std::shared_ptr<CMFDAcceleration>
CMFDAcceleration::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<CMFDAcceleration>("lbs::CMFDAcceleration", params);
}

CMFDAcceleration::CMFDAcceleration(const InputParameters& params)
  : DiscreteOrdinatesKEigenAcceleration(params),
    coarse_mesh_type_(params.GetParamValue<std::string>("coarse_mesh")),
    current_closure_(params.GetParamValue<std::string>("current_closure")),
    automatic_closure_(current_closure_ == "auto"),
    aggregation_size_(params.GetParamValue<int>("aggregation_size")),
    group_aggregation_size_(
      static_cast<unsigned int>(params.GetParamValue<int>("group_aggregation_size"))),
    relaxation_(params.GetParamValue<double>("relaxation")),
    correction_max_attempts_(params.GetParamValue<unsigned int>("correction_max_attempts")),
    correction_min_damping_(params.GetParamValue<double>("correction_min_damping")),
    negative_flux_tolerance_(params.GetParamValue<double>("negative_flux_tolerance")),
    inactive_iterations_(params.GetParamValue<unsigned int>("inactive_iterations")),
    update_wgs_max_its_(params.GetParamValue<unsigned int>("update_wgs_max_its")),
    update_wgs_abs_tol_(params.GetParamValue<double>("update_wgs_abs_tol")),
    balance_residual_tolerance_(params.GetParamValue<double>("balance_residual_tolerance")),
    coarse_solver_policy_(params.GetParamValue<std::string>("coarse_solver_policy")),
    coarse_direct_solve_threshold_(
      static_cast<std::size_t>(params.GetParamValue<int>("coarse_direct_solve_threshold"))),
    num_groups_(do_problem_.GetNumGroups()),
    num_coarse_groups_((num_groups_ + group_aggregation_size_ - 1) / group_aggregation_size_),
    current_closure_blend_(current_closure_ == "partial" ? 1.0 : 0.0)
{
  active_current_closure_ = automatic_closure_ ? "net" : current_closure_;
}

CMFDAcceleration::~CMFDAcceleration()
{
  if (ksp_)
    OpenSnPETScCall(KSPDestroy(&ksp_));
  if (rhs_)
    OpenSnPETScCall(VecDestroy(&rhs_));
  if (A_)
    OpenSnPETScCall(MatDestroy(&A_));
}

void
CMFDAcceleration::Initialize()
{
  do_problem_.ConfigureOutflowStorage(true);
  active_current_closure_ = automatic_closure_ ? "net" : current_closure_;
  current_closure_blend_ = current_closure_ == "partial" ? 1.0 : 0.0;

  if (coarse_mesh_type_ == "identity")
    coarse_mesh_ = CMFDCoarseMesh::BuildIdentity(*do_problem_.GetGrid());
  else if (coarse_mesh_type_ == "local_aggregation")
    coarse_mesh_ = CMFDCoarseMesh::BuildLocalAggregation(*do_problem_.GetGrid(), aggregation_size_);
  else if (coarse_mesh_type_ == "global_aggregation")
    coarse_mesh_ =
      CMFDCoarseMesh::BuildGlobalAggregation(*do_problem_.GetGrid(), aggregation_size_);

  ConfigureTransportSolve(update_wgs_max_its_, update_wgs_abs_tol_);
  log.Log() << no_wrap << program_timer.GetTimeString() << " CMFD configured "
            << do_problem_.GetNumGroupsets()
            << " groupset WGS solvers with max iterations = " << update_wgs_max_its_
            << ", absolute tolerance = " << update_wgs_abs_tol_ << ".";
  InitializeLinearSystem();
  AssembleOperator();

  CoarseMeshDiagnostics diagnostics;
  diagnostics.local_coarse_cells = coarse_mesh_.NumLocalCells();

  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    diagnostics.local_fine_cells += coarse_cell.fine_cell_ids.size();
    diagnostics.local_faces_per_coarse_cell_sum += static_cast<double>(coarse_cell.faces.size());
    if (coarse_cell.fine_cell_ids.size() < static_cast<std::size_t>(aggregation_size_))
      ++diagnostics.local_undersized_coarse_cells;
    diagnostics.local_max_fine_cells_per_coarse_cell =
      std::max(diagnostics.local_max_fine_cells_per_coarse_cell, coarse_cell.fine_cell_ids.size());
    diagnostics.local_max_faces_per_coarse_cell =
      std::max(diagnostics.local_max_faces_per_coarse_cell, coarse_cell.faces.size());
  }

  opensn::mpi_comm.all_reduce(
    diagnostics.local_fine_cells, diagnostics.global_fine_cells, mpi::op::sum<std::size_t>());
  opensn::mpi_comm.all_reduce(
    diagnostics.local_coarse_cells, diagnostics.global_coarse_cells, mpi::op::sum<std::size_t>());
  opensn::mpi_comm.all_reduce(diagnostics.local_max_fine_cells_per_coarse_cell,
                              diagnostics.global_max_fine_cells_per_coarse_cell,
                              mpi::op::max<std::size_t>());
  opensn::mpi_comm.all_reduce(diagnostics.local_max_faces_per_coarse_cell,
                              diagnostics.global_max_faces_per_coarse_cell,
                              mpi::op::max<std::size_t>());
  opensn::mpi_comm.all_reduce(diagnostics.local_undersized_coarse_cells,
                              diagnostics.global_undersized_coarse_cells,
                              mpi::op::sum<std::size_t>());
  opensn::mpi_comm.all_reduce(diagnostics.local_faces_per_coarse_cell_sum,
                              diagnostics.global_faces_per_coarse_cell_sum,
                              mpi::op::sum<double>());

  const double aggregation_ratio = diagnostics.global_coarse_cells > 0
                                     ? static_cast<double>(diagnostics.global_fine_cells) /
                                         static_cast<double>(diagnostics.global_coarse_cells)
                                     : 0.0;
  const double average_faces_per_coarse_cell =
    diagnostics.global_coarse_cells > 0 ? diagnostics.global_faces_per_coarse_cell_sum /
                                            static_cast<double>(diagnostics.global_coarse_cells)
                                        : 0.0;

  log.Log() << no_wrap << program_timer.GetTimeString() << " Initialized CMFD coarse mesh with "
            << diagnostics.global_coarse_cells << " global coarse cells from "
            << diagnostics.global_fine_cells << " fine cells"
            << " (avg fine/coarse=" << aggregation_ratio << ", coarse groups=" << num_coarse_groups_
            << ", global unknowns=" << GlobalUnknownCount() << ").";
  if (verbose_)
  {
    const auto LogMetric = [](const std::string& metric, const auto value)
    { LogCMFDAcceleratorMetric("coarse_mesh", -1, metric, MetricValue(value)); };
    LogMetric("global_fine_cells", diagnostics.global_fine_cells);
    LogMetric("global_coarse_cells", diagnostics.global_coarse_cells);
    LogMetric("group_aggregation_size", group_aggregation_size_);
    LogMetric("energy_coarse_groups", num_coarse_groups_);
    LogMetric("aggregation_ratio", aggregation_ratio);
    LogMetric("global_unknowns", GlobalUnknownCount());
    LogMetric("max_fine_cells_per_coarse_cell", diagnostics.global_max_fine_cells_per_coarse_cell);
    LogMetric("max_faces_per_coarse_cell", diagnostics.global_max_faces_per_coarse_cell);
    LogMetric("local_coarse_cells", diagnostics.local_coarse_cells);
    LogMetric("undersized_coarse_cells", diagnostics.global_undersized_coarse_cells);
    LogMetric("average_faces_per_coarse_cell", average_faces_per_coarse_cell);
  }
}

void
CMFDAcceleration::ProbeAutomaticClosure(const double k_eff,
                                        const std::vector<double>& coarse_phi,
                                        const std::vector<double>& coarse_source_phi)
{
  if (not automatic_closure_ or automatic_closure_probed_)
    return;
  if (automatic_closure_probe_iterations_ >= 8)
  {
    automatic_closure_probed_ = true;
    return;
  }

  const auto previous_closure = active_current_closure_;
  const double previous_blend = current_closure_blend_;
  const auto SolveCoarseK = [this, k_eff, &coarse_phi]()
  {
    auto coarse_phi_m = coarse_phi;
    auto coarse_phi_candidate = coarse_phi;
    double coarse_k_eff = k_eff;
    double coarse_production_old = ComputeCoarseFissionProduction(coarse_phi_m);
    OpenSnLogicalErrorIf(coarse_production_old == 0.0,
                         "CMFDAcceleration cannot probe with zero coarse production.");

    const int max_probe_iterations = std::min<int>(pi_max_its_, 8);
    for (int m = 0; m < max_probe_iterations; ++m)
    {
      coarse_phi_candidate = SolveCoarseSystem(coarse_k_eff, coarse_phi_m);
      const double coarse_production_new = ComputeCoarseFissionProduction(coarse_phi_candidate);
      OpenSnLogicalErrorIf(coarse_production_new == 0.0,
                           "CMFDAcceleration cannot update probed k with zero production.");
      const double coarse_k_eff_new = coarse_production_new / coarse_production_old * coarse_k_eff;
      const double k_change = std::fabs(coarse_k_eff_new / coarse_k_eff - 1.0);
      coarse_k_eff = coarse_k_eff_new;
      coarse_phi_m = coarse_phi_candidate;
      coarse_production_old = coarse_production_new;
      if (k_change < pi_k_tol_)
        break;
    }
    return coarse_k_eff;
  };

  active_current_closure_ = "net";
  current_closure_blend_ = 0.0;
  AssembleOperator();
  const auto net_residual = ComputeBalanceResidual(k_eff, coarse_phi, coarse_source_phi);
  const double net_coarse_k_eff = SolveCoarseK();

  active_current_closure_ = "partial";
  current_closure_blend_ = 1.0;
  AssembleOperator();
  const auto partial_residual = ComputeBalanceResidual(k_eff, coarse_phi, coarse_source_phi);
  const double partial_coarse_k_eff = SolveCoarseK();

  const double net_l2 = net_residual.relative_l2;
  const double partial_l2 = partial_residual.relative_l2;
  const double net_k_change =
    std::fabs(net_coarse_k_eff - k_eff) / std::max(std::fabs(k_eff), 1.0e-30);
  const double partial_k_change =
    std::fabs(partial_coarse_k_eff - k_eff) / std::max(std::fabs(k_eff), 1.0e-30);
  const double materiality = std::max(1.0e-3, 5.0 * balance_residual_tolerance_);
  const double residual_ratio = net_l2 > materiality ? net_l2 / std::max(partial_l2, 1.0e-30) : 1.0;
  const double k_ratio =
    net_k_change > 1.0e-2 ? net_k_change / std::max(partial_k_change, 1.0e-30) : 1.0;
  const double quality_ratio = std::max(residual_ratio, k_ratio);
  const double selected_blend =
    quality_ratio <= 1.1 ? 0.0 : std::clamp((quality_ratio - 1.1) / 0.9, 0.0, 1.0);

  active_current_closure_ =
    selected_blend >= 1.0 ? "partial" : (selected_blend > 0.0 ? "blend" : previous_closure);
  current_closure_blend_ = selected_blend > 0.0 ? selected_blend : previous_blend;
  automatic_closure_probe_iterations_ += 1;
  automatic_closure_probed_ = selected_blend > 0.0;

  if (selected_blend > 0.0)
  {
    log.Log0() << no_wrap << program_timer.GetTimeString() << " CMFD auto closure selected "
               << active_current_closure_ << ": partial fraction = " << selected_blend
               << ", net residual = " << net_l2 << ", partial residual = " << partial_l2
               << ", net coarse k = " << net_coarse_k_eff
               << ", partial coarse k = " << partial_coarse_k_eff << ".";
  }
  else if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
    log.Log() << no_wrap << program_timer.GetTimeString() << " CMFD auto closure probe kept net"
              << ": net residual = " << net_l2 << ", partial residual = " << partial_l2
              << ", net coarse k = " << net_coarse_k_eff
              << ", partial coarse k = " << partial_coarse_k_eff
              << ", selected blend = " << selected_blend
              << ", probe = " << automatic_closure_probe_iterations_ << ".";
}

void
CMFDAcceleration::UpdateAutomaticClosure(const BalanceResidual operator_residual,
                                         const BalanceResidual transport_current_residual,
                                         const double transport_k_eff,
                                         const double coarse_k_eff,
                                         const double applied_damping,
                                         const bool skipped_correction)
{
  if (not automatic_closure_ or active_current_closure_ != "net")
    return;

  const double residual_ratio =
    operator_residual.relative_l2 / std::max(transport_current_residual.relative_l2, 1.0e-30);
  const double coarse_k_change =
    std::fabs(coarse_k_eff - transport_k_eff) / std::max(std::fabs(transport_k_eff), 1.0e-30);
  const double expected_damping = std::max(relaxation_, 1.0e-30);
  const bool strongly_damped = applied_damping < 0.25 * expected_damping;
  const bool startup_iteration = outer_iteration_ < 8;
  const bool bad_balance =
    startup_iteration and residual_ratio > 5.0 and
    operator_residual.relative_l2 > std::max(1.0e-3, 5.0 * balance_residual_tolerance_);
  const bool bad_coarse_k = startup_iteration and coarse_k_change > 0.02;
  const bool bad_correction = skipped_correction or strongly_damped;

  if (bad_balance and (bad_coarse_k or bad_correction or residual_ratio > 25.0))
  {
    active_current_closure_ = "partial";

    log.Log0() << no_wrap << program_timer.GetTimeString()
               << " CMFD auto closure switched from net to partial"
               << ": operator/current residual ratio = " << residual_ratio
               << ", coarse_k_change = " << coarse_k_change
               << ", applied_damping = " << applied_damping << ".";
  }
}

void
CMFDAcceleration::PreExecute()
{
  coarse_phi_old_fine_.clear();
  coarse_phi_new_fine_.clear();
  coarse_phi_old_.clear();
  coarse_phi_new_.clear();
  coarse_phi_cache_.clear();
  face_current_cache_.clear();
  outer_iteration_ = 0;
  active_current_closure_ = automatic_closure_ ? "net" : current_closure_;
  current_closure_blend_ = current_closure_ == "partial" ? 1.0 : 0.0;
  automatic_closure_probed_ = false;
  automatic_closure_probe_iterations_ = 0;
  ConfigureTransportSolve(update_wgs_max_its_, update_wgs_abs_tol_);
  last_restrict_old_time_ = 0.0;
  transport_start_time_ = 0.0;
  old_fission_production_ = 0.0;
  consecutive_skipped_corrections_ = 0;
}

void
CMFDAcceleration::PrePowerIteration()
{
  const double t0 = NowSeconds();
  old_fission_production_ = ComputeFissionProduction(do_problem_, phi_old_local_);
  OpenSnLogicalErrorIf(old_fission_production_ == 0.0,
                       "CMFDAcceleration cannot update k with zero old production.");
  coarse_phi_old_fine_ =
    CMFDRestrictScalarFlux(do_problem_, first_group_, num_groups_, coarse_mesh_, phi_old_local_);
  coarse_phi_old_ = CMFDRestrictScalarFlux(
    do_problem_, first_group_, num_groups_, group_aggregation_size_, coarse_mesh_, phi_old_local_);
  last_restrict_old_time_ = NowSeconds() - t0;
  transport_start_time_ = NowSeconds();
}

double
CMFDAcceleration::PostPowerIteration()
{
  const double post_start_time = NowSeconds();
  const double transport_time =
    transport_start_time_ > 0.0 ? post_start_time - transport_start_time_ : 0.0;

  double t0 = NowSeconds();
  coarse_phi_new_fine_ =
    CMFDRestrictScalarFlux(do_problem_, first_group_, num_groups_, coarse_mesh_, phi_new_local_);
  coarse_phi_new_ = CMFDRestrictScalarFlux(
    do_problem_, first_group_, num_groups_, group_aggregation_size_, coarse_mesh_, phi_new_local_);
  const double restrict_new_time = NowSeconds() - t0;

  t0 = NowSeconds();
  const double F_transport_new = ComputeFissionProduction(do_problem_, phi_new_local_);
  OpenSnLogicalErrorIf(old_fission_production_ == 0.0,
                       "CMFDAcceleration cannot update k with zero old production.");
  const double raw_transport_k_eff =
    F_transport_new / old_fission_production_ * solver_->GetEigenvalue();
  const double cmfd_reference_fission_production =
    ComputeFissionProduction(do_problem_, phi_old_local_);
  OpenSnLogicalErrorIf(cmfd_reference_fission_production == 0.0,
                       "CMFDAcceleration cannot update k with zero CMFD reference production.");
  const double transport_k_eff =
    F_transport_new / cmfd_reference_fission_production * solver_->GetEigenvalue();
  const double fission_time = NowSeconds() - t0;

  t0 = NowSeconds();
  BuildCoarseFluxCache();
  const double coarse_flux_cache_time = NowSeconds() - t0;

  t0 = NowSeconds();
  BuildFaceCurrentCache();
  const double face_current_cache_time = NowSeconds() - t0;

  if (outer_iteration_ < inactive_iterations_)
  {
    for (const auto& groupset : groupsets_)
      LBSVecOps::GSScopedCopyPrimarySTLvectors(
        do_problem_, groupset, phi_new_local_, phi_old_local_);
    last_update_allows_convergence_ = false;
    last_balance_residual_known_ = false;
    last_correction_skipped_ = false;
    last_correction_skip_reason_.clear();
    if (do_problem_.GetOptions().verbose_outer_iterations)
    {
      std::ostringstream iter_info;
      iter_info << program_timer.GetTimeString() << " CMFD iteration = " << outer_iteration_;
      AppendNumericField(iter_info, "k_eff", transport_k_eff, Fixed(7));
      iter_info << ", status = inactive";
      AppendNumericField(iter_info, "inactive_iterations", inactive_iterations_);
      log.Log() << no_wrap << iter_info.str();
    }

    ++outer_iteration_;
    return transport_k_eff;
  }

  ProbeAutomaticClosure(transport_k_eff, coarse_phi_new_, coarse_phi_old_);

  t0 = NowSeconds();
  AssembleOperator();
  const double assemble_time = NowSeconds() - t0;

  t0 = NowSeconds();
  const auto transport_balance_residual =
    ComputeBalanceResidual(transport_k_eff, coarse_phi_new_, coarse_phi_old_);
  const auto transport_current_balance_residual =
    ComputeTransportCurrentBalanceResidual(transport_k_eff, coarse_phi_new_, coarse_phi_old_);
  const bool transport_balance_allows_convergence =
    transport_current_balance_residual.relative_l2 <= balance_residual_tolerance_;
  last_balance_residual_known_ = true;
  last_transport_current_balance_relative_l2_ = transport_current_balance_residual.relative_l2;
  const double residual_time = NowSeconds() - t0;
  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
    log.Log() << no_wrap << program_timer.GetTimeString()
              << " CMFD operator balance residual: max_abs = " << transport_balance_residual.max_abs
              << ", relative_l2 = " << transport_balance_residual.relative_l2;
  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
    log.Log() << no_wrap << program_timer.GetTimeString()
              << " CMFD transport-current balance residual: max_abs = "
              << transport_current_balance_residual.max_abs
              << ", relative_l2 = " << transport_current_balance_residual.relative_l2;

  auto coarse_phi_m = coarse_phi_new_;
  auto coarse_phi = coarse_phi_new_;
  double k_eff = transport_k_eff;
  double coarse_production_old = ComputeCoarseFissionProduction(coarse_phi_m);
  OpenSnLogicalErrorIf(coarse_production_old == 0.0,
                       "CMFDAcceleration cannot iterate with zero coarse production.");
  int coarse_pi_iterations = 0;
  int coarse_ksp_iterations = 0;
  bool coarse_linear_solves_converged = true;
  t0 = NowSeconds();
  for (int m = 0; m < pi_max_its_; ++m)
  {
    coarse_phi = SolveCoarseSystem(k_eff, coarse_phi_m);
    coarse_ksp_iterations += last_ksp_iterations_;
    coarse_linear_solves_converged = coarse_linear_solves_converged and last_ksp_converged_;
    coarse_pi_iterations = m + 1;
    const double coarse_production_new = ComputeCoarseFissionProduction(coarse_phi);
    OpenSnLogicalErrorIf(coarse_production_new == 0.0,
                         "CMFDAcceleration cannot update k with zero coarse production.");
    const double k_eff_new = coarse_production_new / coarse_production_old * k_eff;
    const double k_change = std::fabs(k_eff_new / k_eff - 1.0);

    k_eff = k_eff_new;
    coarse_phi_m = coarse_phi;
    coarse_production_old = coarse_production_new;

    if (k_change < pi_k_tol_)
      break;
  }
  const double coarse_pi_time = NowSeconds() - t0;
  const double coarse_phi_change = RelativeL2Difference(coarse_phi, coarse_phi_new_);
  const double coarse_k_eff = k_eff;

  t0 = NowSeconds();
  double applied_damping = 0.0;
  unsigned int correction_attempts = 0;
  FluxUpdateDiagnostics correction_diagnostics;
  bool skipped_correction = false;
  std::string correction_skip_reason;
  if (not coarse_linear_solves_converged)
  {
    skipped_correction = true;
    correction_skip_reason = "coarse_linear_solve_not_converged";
    correction_diagnostics = AnalyzeFluxUpdate(phi_new_local_, transport_k_eff);
    k_eff = transport_k_eff;
    const bool can_switch_to_direct = coarse_solver_policy_ == "auto" and
                                      not using_direct_coarse_solver_ and
                                      GlobalUnknownCount() <= coarse_direct_solve_threshold_;
    if (can_switch_to_direct)
    {
      ConfigureCoarseLinearSolver(true);
      using_direct_coarse_solver_ = true;
    }
    if (verbose_ or outer_iteration_ < 5 or outer_iteration_ % 25 == 0)
    {
      std::ostringstream warning;
      warning << "CMFD correction skipped because one or more coarse linear solves did not "
                 "converge.";
      if (can_switch_to_direct)
        warning << " Switching auto coarse solver policy to direct for subsequent corrections.";
      log.Log0Warning() << no_wrap << warning.str();
    }
  }
  else
    k_eff = ApplyFluxCorrectionWithDamping(phi_new_local_,
                                           coarse_phi,
                                           cmfd_reference_fission_production,
                                           transport_k_eff,
                                           applied_damping,
                                           correction_attempts,
                                           correction_diagnostics,
                                           skipped_correction,
                                           correction_skip_reason);
  if (skipped_correction)
  {
    ++consecutive_skipped_corrections_;
    if (consecutive_skipped_corrections_ >= 2)
      k_eff = raw_transport_k_eff;
  }
  else
    consecutive_skipped_corrections_ = 0;
  last_correction_skipped_ = skipped_correction;
  last_correction_skip_reason_ = skipped_correction ? correction_skip_reason : std::string();
  last_update_allows_convergence_ = transport_balance_allows_convergence and not skipped_correction;
  const double accelerated_k_change = std::fabs(k_eff - solver_->GetEigenvalue()) / k_eff;
  UpdateAutomaticClosure(transport_balance_residual,
                         transport_current_balance_residual,
                         transport_k_eff,
                         coarse_k_eff,
                         applied_damping,
                         skipped_correction);
  for (const auto& groupset : groupsets_)
    LBSVecOps::GSScopedCopyPrimarySTLvectors(do_problem_, groupset, phi_new_local_, phi_old_local_);
  const double update_time = NowSeconds() - t0;
  const double post_time = NowSeconds() - post_start_time;

  if (do_problem_.GetOptions().verbose_outer_iterations)
  {
    std::ostringstream iter_info;
    iter_info << program_timer.GetTimeString() << " CMFD iteration = " << outer_iteration_;
    AppendNumericField(iter_info, "k_eff", k_eff, Fixed(7));
    iter_info << ", status = " << (skipped_correction ? "skipped" : "applied");
    log.Log() << no_wrap << iter_info.str();
  }

  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
  {
    const auto LogMetric = [this](const std::string& metric, const auto value)
    { LogCMFDAcceleratorMetric("correction", outer_iteration_, metric, MetricValue(value)); };
    LogMetric("relaxation", relaxation_);
    LogMetric("damping", applied_damping);
    LogMetric("attempts", correction_attempts);
    LogMetric("skipped", skipped_correction ? 1 : 0);
    LogMetric("min_scalar_flux", correction_diagnostics.min_scalar_flux);
    LogMetric("nonfinite_flux", correction_diagnostics.has_nonfinite_flux ? 1 : 0);
    LogMetric("invalid_k", correction_diagnostics.has_invalid_k ? 1 : 0);
    LogMetric("transport_k_eff", transport_k_eff);
    LogMetric("coarse_k_eff", coarse_k_eff);
    LogMetric("applied_k_eff", k_eff);
    LogMetric("coarse_phi_change", coarse_phi_change);
  }

  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
  {
    const auto LogMetric = [this](const std::string& metric, const auto value)
    { LogCMFDAcceleratorMetric("timing", outer_iteration_, metric, MetricValue(value)); };
    LogMetric("transport_time", MaxAcrossRanks(transport_time));
    LogMetric("restrict_old_time", MaxAcrossRanks(last_restrict_old_time_));
    LogMetric("restrict_new_time", MaxAcrossRanks(restrict_new_time));
    LogMetric("coarse_flux_cache_time", MaxAcrossRanks(coarse_flux_cache_time));
    LogMetric("face_current_cache_time", MaxAcrossRanks(face_current_cache_time));
    LogMetric("assemble_time", MaxAcrossRanks(assemble_time));
    LogMetric("residual_time", MaxAcrossRanks(residual_time));
    LogMetric("fission_time", MaxAcrossRanks(fission_time));
    LogMetric("coarse_pi_time", MaxAcrossRanks(coarse_pi_time));
    LogMetric("update_time", MaxAcrossRanks(update_time));
    LogMetric("post_time", MaxAcrossRanks(post_time));
    LogMetric("coarse_pi_iterations", coarse_pi_iterations);
    LogMetric("coarse_ksp_iterations", coarse_ksp_iterations);
    LogMetric("coarse_linear_solves_converged", coarse_linear_solves_converged ? 1 : 0);
    LogMetric("operator_balance_residual", transport_balance_residual.relative_l2);
    LogMetric("transport_current_balance_residual", transport_current_balance_residual.relative_l2);
    LogMetric("balance_residual_tolerance", balance_residual_tolerance_);
    LogMetric("balance_allows_convergence", transport_balance_allows_convergence ? 1 : 0);
    LogMetric("accelerated_k_change", accelerated_k_change);
  }

  ++outer_iteration_;
  return k_eff;
}

std::string
CMFDAcceleration::GetPowerIterationConvergenceInfo() const
{
  if (not last_balance_residual_known_)
    return {};

  std::ostringstream out;
  out << "CMFD convergence check = "
      << (last_update_allows_convergence_ ? "satisfied" : "not_satisfied")
      << ", transport-current balance residual = " << std::scientific << std::setprecision(6)
      << last_transport_current_balance_relative_l2_
      << " (balance_residual_tolerance = " << balance_residual_tolerance_ << ")";
  out << std::defaultfloat;
  if (last_correction_skipped_)
  {
    out << ", correction = skipped";
    if (not last_correction_skip_reason_.empty())
      out << " (" << last_correction_skip_reason_ << ")";
  }
  return out.str();
}

std::size_t
CMFDAcceleration::LocalUnknownCount() const
{
  return coarse_mesh_.NumLocalCells() * num_coarse_groups_;
}

std::size_t
CMFDAcceleration::GlobalUnknownCount() const
{
  return coarse_mesh_.NumGlobalCells() * num_coarse_groups_;
}

PetscInt
CMFDAcceleration::MapDOF(const uint64_t coarse_cell_global_id,
                         const unsigned int coarse_group) const
{
  return static_cast<PetscInt>(coarse_cell_global_id * num_coarse_groups_ + coarse_group);
}

unsigned int
CMFDAcceleration::FineGroupBegin(const unsigned int coarse_group) const
{
  return coarse_group * group_aggregation_size_;
}

unsigned int
CMFDAcceleration::FineGroupEnd(const unsigned int coarse_group) const
{
  return std::min(num_groups_, (coarse_group + 1) * group_aggregation_size_);
}

double
CMFDAcceleration::CondensedRemovalXS(const MultiGroupXS& xs,
                                     const std::vector<double>& fine_coarse_phi,
                                     const CMFDCoarseCell& coarse_cell,
                                     const unsigned int coarse_group) const
{
  const auto& sigma_t = xs.GetSigmaTotal();
  const auto& scatter = xs.GetTransferMatrix(0);
  const bool has_flux = not fine_coarse_phi.empty();
  const auto fine_offset = coarse_cell.local_id * num_groups_;
  double denominator = 0.0;
  for (unsigned int g = FineGroupBegin(coarse_group); g < FineGroupEnd(coarse_group); ++g)
    denominator += has_flux ? fine_coarse_phi[fine_offset + g] : 1.0;

  double numerator = 0.0;
  for (unsigned int gp = FineGroupBegin(coarse_group); gp < FineGroupEnd(coarse_group); ++gp)
    numerator += sigma_t[gp] * (has_flux ? fine_coarse_phi[fine_offset + gp] : 1.0);

  for (unsigned int g = FineGroupBegin(coarse_group); g < FineGroupEnd(coarse_group); ++g)
    for (const auto& [_, gp, sigma_s] : scatter.Row(g))
      if (gp >= FineGroupBegin(coarse_group) and gp < FineGroupEnd(coarse_group))
        numerator -= sigma_s * (has_flux ? fine_coarse_phi[fine_offset + gp] : 1.0);

  if (std::fabs(denominator) <= 1.0e-30)
  {
    numerator = 0.0;
    denominator = 0.0;
    for (unsigned int gp = FineGroupBegin(coarse_group); gp < FineGroupEnd(coarse_group); ++gp)
    {
      numerator += sigma_t[gp];
      denominator += 1.0;
    }
    for (unsigned int g = FineGroupBegin(coarse_group); g < FineGroupEnd(coarse_group); ++g)
      for (const auto& [_, gp, sigma_s] : scatter.Row(g))
        if (gp >= FineGroupBegin(coarse_group) and gp < FineGroupEnd(coarse_group))
          numerator -= sigma_s;
  }
  return numerator / denominator;
}

double
CMFDAcceleration::CondensedDiffusionCoefficient(const MultiGroupXS& xs,
                                                const std::vector<double>& fine_coarse_phi,
                                                const CMFDCoarseCell& coarse_cell,
                                                const unsigned int coarse_group) const
{
  const auto& diffusion_coeff = xs.GetDiffusionCoefficient();
  const bool has_flux = not fine_coarse_phi.empty();
  const auto fine_offset = coarse_cell.local_id * num_groups_;
  double numerator = 0.0;
  double denominator = 0.0;
  for (unsigned int g = FineGroupBegin(coarse_group); g < FineGroupEnd(coarse_group); ++g)
  {
    const double phi = has_flux ? fine_coarse_phi[fine_offset + g] : 1.0;
    numerator += diffusion_coeff[g] * phi;
    denominator += phi;
  }
  if (std::fabs(denominator) <= 1.0e-30)
  {
    numerator = 0.0;
    denominator = 0.0;
    for (unsigned int g = FineGroupBegin(coarse_group); g < FineGroupEnd(coarse_group); ++g)
    {
      numerator += diffusion_coeff[g];
      denominator += 1.0;
    }
  }
  return numerator / denominator;
}

double
CMFDAcceleration::CondensedScatterXS(const MultiGroupXS& xs,
                                     const std::vector<double>& fine_coarse_phi,
                                     const CMFDCoarseCell& coarse_cell,
                                     const unsigned int row_coarse_group,
                                     const unsigned int col_coarse_group) const
{
  const auto& scatter = xs.GetTransferMatrix(0);
  const bool has_flux = not fine_coarse_phi.empty();
  const auto fine_offset = coarse_cell.local_id * num_groups_;
  double numerator = 0.0;
  double denominator = 0.0;
  for (unsigned int gp = FineGroupBegin(col_coarse_group); gp < FineGroupEnd(col_coarse_group);
       ++gp)
    denominator += has_flux ? fine_coarse_phi[fine_offset + gp] : 1.0;

  for (unsigned int g = FineGroupBegin(row_coarse_group); g < FineGroupEnd(row_coarse_group); ++g)
    for (const auto& [_, gp, sigma_s] : scatter.Row(g))
      if (gp >= FineGroupBegin(col_coarse_group) and gp < FineGroupEnd(col_coarse_group))
        numerator += sigma_s * (has_flux ? fine_coarse_phi[fine_offset + gp] : 1.0);

  if (std::fabs(denominator) <= 1.0e-30)
  {
    numerator = 0.0;
    denominator = 0.0;
    for (unsigned int gp = FineGroupBegin(col_coarse_group); gp < FineGroupEnd(col_coarse_group);
         ++gp)
      denominator += 1.0;
    for (unsigned int g = FineGroupBegin(row_coarse_group); g < FineGroupEnd(row_coarse_group); ++g)
      for (const auto& [_, gp, sigma_s] : scatter.Row(g))
        if (gp >= FineGroupBegin(col_coarse_group) and gp < FineGroupEnd(col_coarse_group))
          numerator += sigma_s;
  }
  return numerator / denominator;
}

double
CMFDAcceleration::CondensedProductionXS(const MultiGroupXS& xs,
                                        const std::vector<double>& fine_coarse_phi,
                                        const CMFDCoarseCell& coarse_cell,
                                        const unsigned int row_coarse_group,
                                        const unsigned int col_coarse_group) const
{
  const auto& F = xs.GetProductionMatrix();
  const bool include_delayed_production = UseDelayedNeutronProduction(do_problem_);
  const bool has_flux = not fine_coarse_phi.empty();
  const auto fine_offset = coarse_cell.local_id * num_groups_;
  double numerator = 0.0;
  double denominator = 0.0;
  for (unsigned int gp = FineGroupBegin(col_coarse_group); gp < FineGroupEnd(col_coarse_group);
       ++gp)
  {
    const double phi = has_flux ? fine_coarse_phi[fine_offset + gp] : 1.0;
    denominator += phi;
    for (unsigned int g = FineGroupBegin(row_coarse_group); g < FineGroupEnd(row_coarse_group); ++g)
    {
      numerator += F[g][gp] * phi;
      if (include_delayed_production)
        numerator += ComputeDelayedFissionProduction(xs, g, gp) * phi;
    }
  }
  if (std::fabs(denominator) <= 1.0e-30)
  {
    numerator = 0.0;
    denominator = 0.0;
    for (unsigned int gp = FineGroupBegin(col_coarse_group); gp < FineGroupEnd(col_coarse_group);
         ++gp)
    {
      denominator += 1.0;
      for (unsigned int g = FineGroupBegin(row_coarse_group); g < FineGroupEnd(row_coarse_group);
           ++g)
      {
        numerator += F[g][gp];
        if (include_delayed_production)
          numerator += ComputeDelayedFissionProduction(xs, g, gp);
      }
    }
  }
  return numerator / denominator;
}

void
CMFDAcceleration::InitializeLinearSystem()
{
  const auto local_unknowns = static_cast<PetscInt>(LocalUnknownCount());
  const auto global_unknowns = static_cast<PetscInt>(GlobalUnknownCount());

  OpenSnPETScCall(MatCreate(opensn::mpi_comm, &A_));
  OpenSnPETScCall(MatSetType(A_, MATMPIAIJ));
  OpenSnPETScCall(
    MatSetSizes(A_, local_unknowns, local_unknowns, global_unknowns, global_unknowns));
  OpenSnPETScCall(MatMPIAIJSetPreallocation(A_, 16, nullptr, 16, nullptr));
  OpenSnPETScCall(MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
  OpenSnPETScCall(MatSetOption(A_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
  OpenSnPETScCall(MatSetUp(A_));

  OpenSnPETScCall(VecCreate(opensn::mpi_comm, &rhs_));
  OpenSnPETScCall(VecSetType(rhs_, VECMPI));
  OpenSnPETScCall(VecSetSizes(rhs_, local_unknowns, global_unknowns));

  OpenSnPETScCall(KSPCreate(opensn::mpi_comm, &ksp_));
  OpenSnPETScCall(KSPSetOperators(ksp_, A_, A_));
  OpenSnPETScCall(KSPSetOptionsPrefix(ksp_, GetName().c_str()));
  const bool petsc_options_policy = coarse_solver_policy_ == "petsc_options";
  using_direct_coarse_solver_ =
    coarse_solver_policy_ == "direct" or
    (coarse_solver_policy_ == "auto" and GlobalUnknownCount() <= coarse_direct_solve_threshold_);
  if (not petsc_options_policy)
    ConfigureCoarseLinearSolver(using_direct_coarse_solver_);
  OpenSnPETScCall(KSPSetTolerances(ksp_, l_abs_tol_, PETSC_DEFAULT, PETSC_DEFAULT, max_iters_));
  if ((petsc_options_policy or petsc_options_ != "ssss") and not petsc_options_.empty())
    OpenSnPETScCall(PetscOptionsInsertString(nullptr, petsc_options_.c_str()));
  OpenSnPETScCall(KSPSetFromOptions(ksp_));
}

void
CMFDAcceleration::ConfigureCoarseLinearSolver(const bool direct)
{
  PC pc = nullptr;
  OpenSnPETScCall(KSPGetPC(ksp_, &pc));
  if (direct)
  {
    OpenSnPETScCall(KSPSetType(ksp_, KSPPREONLY));
    OpenSnPETScCall(PCSetType(pc, PCLU));
  }
  else
  {
    OpenSnPETScCall(KSPSetType(ksp_, KSPGMRES));
    OpenSnPETScCall(PCSetType(pc, PCJACOBI));
  }
}

void
CMFDAcceleration::ConfigureTransportSolve(const unsigned int max_iterations,
                                          const double residual_tolerance)
{
  for (std::size_t gsid = 0; gsid < do_problem_.GetNumGroupsets(); ++gsid)
  {
    auto& groupset = do_problem_.GetGroupset(gsid);
    groupset.max_iterations = max_iterations;
    groupset.residual_tolerance = residual_tolerance;

    auto wgs_solver = do_problem_.GetWGSSolver(groupset.id);
    if (auto petsc_solver = std::dynamic_pointer_cast<PETScLinearSolver>(wgs_solver))
    {
      auto& tolerance_options = petsc_solver->GetToleranceOptions();
      tolerance_options.maximum_iterations = static_cast<PetscInt>(max_iterations);
      tolerance_options.residual_absolute = residual_tolerance;
      petsc_solver->ApplyToleranceOptions();
    }
    else
    {
      OpenSnLogicalErrorIf(
        not std::dynamic_pointer_cast<ClassicRichardson>(wgs_solver),
        "CMFD transport updates require a PETSc or classic Richardson WGS solver.");
    }
  }
}

void
CMFDAcceleration::BuildCoarseFluxCache()
{
  coarse_phi_cache_.clear();

  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      coarse_phi_cache_[std::make_tuple(coarse_cell.global_id, cg)] = coarse_phi_new_[offset + cg];
  }

  std::map<int, std::set<std::pair<uint64_t, unsigned int>>> pid_request_sets;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    for (const auto& face : coarse_cell.faces)
    {
      if (not face.has_neighbor or face.neighbor_partition_id == opensn::mpi_comm.rank())
        continue;

      auto& requests = pid_request_sets[face.neighbor_partition_id];
      for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
        requests.emplace(face.neighbor_id, cg);
    }
  }

  std::map<int, std::vector<uint64_t>> pid_requests;
  for (const auto& [pid, request_set] : pid_request_sets)
  {
    auto& requests = pid_requests[pid];
    requests.reserve(2 * request_set.size());
    for (const auto& [coarse_cell_gid, g] : request_set)
    {
      requests.push_back(coarse_cell_gid);
      requests.push_back(g);
    }
  }

  const auto received_requests = MapAllToAll(pid_requests);

  std::map<int, std::vector<uint64_t>> pid_response_keys;
  std::map<int, std::vector<double>> pid_response_values;
  for (const auto& [pid, requests] : received_requests)
  {
    OpenSnLogicalErrorIf(requests.size() % 2 != 0, "Invalid CMFD coarse-flux request buffer.");
    auto& response_keys = pid_response_keys[pid];
    auto& response_values = pid_response_values[pid];
    response_keys.reserve(requests.size());
    response_values.reserve(requests.size() / 2);

    for (std::size_t r = 0; r < requests.size(); r += 2)
    {
      const auto coarse_cell_gid = requests[r];
      const auto cg = static_cast<unsigned int>(requests[r + 1]);
      OpenSnInvalidArgumentIf(not coarse_mesh_.HasLocalCoarseCell(coarse_cell_gid),
                              "CMFD coarse-flux request references a nonlocal coarse cell.");

      const auto coarse_local_id = coarse_mesh_.LocalCellFromGlobalID(coarse_cell_gid).local_id;
      const auto value = coarse_phi_new_[coarse_local_id * num_coarse_groups_ + cg];
      response_keys.push_back(coarse_cell_gid);
      response_keys.push_back(cg);
      response_values.push_back(value);
    }
  }

  const auto received_response_keys = MapAllToAll(pid_response_keys);
  const auto received_response_values = MapAllToAll(pid_response_values);
  for (const auto& [pid, keys] : received_response_keys)
  {
    const auto values_it = received_response_values.find(pid);
    OpenSnLogicalErrorIf(values_it == received_response_values.end(),
                         "Missing CMFD coarse-flux response values.");
    const auto& values = values_it->second;
    OpenSnLogicalErrorIf(keys.size() != 2 * values.size(),
                         "Invalid CMFD coarse-flux response buffer.");

    for (std::size_t r = 0; r < values.size(); ++r)
    {
      const auto key_offset = 2 * r;
      const auto coarse_cell_gid = keys[key_offset];
      const auto cg = static_cast<unsigned int>(keys[key_offset + 1]);
      coarse_phi_cache_[std::make_tuple(coarse_cell_gid, cg)] = values[r];
    }
  }
}

void
CMFDAcceleration::BuildFaceCurrentCache()
{
  const auto& grid = *do_problem_.GetGrid();
  const auto& outflow_views = do_problem_.GetCellOutflowViews();

  std::map<int, std::set<std::tuple<uint64_t, std::size_t, uint64_t, unsigned int>>>
    pid_request_sets;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    for (const auto& face : coarse_cell.faces)
    {
      for (const auto& fine_face : face.fine_faces)
      {
        if (fine_face.cell_partition_id != opensn::mpi_comm.rank())
        {
          auto& requests = pid_request_sets[fine_face.cell_partition_id];
          for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
            requests.emplace(fine_face.cell_id,
                             fine_face.face_index,
                             fine_face.neighbor_id.value_or(std::numeric_limits<uint64_t>::max()),
                             cg);
        }

        const auto neighbor_id = fine_face.neighbor_id;
        if (neighbor_id.has_value() and fine_face.neighbor_partition_id != opensn::mpi_comm.rank())
        {
          auto& requests = pid_request_sets[fine_face.neighbor_partition_id];
          for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
            requests.emplace(
              *neighbor_id, std::numeric_limits<std::size_t>::max(), fine_face.cell_id, cg);
        }
      }
    }
  }

  std::map<int, std::vector<uint64_t>> pid_requests;
  for (const auto& [pid, request_set] : pid_request_sets)
  {
    auto& requests = pid_requests[pid];
    requests.reserve(4 * request_set.size());
    for (const auto& [owner_cell_gid, owner_face_index, neighbor_cell_gid, cg] : request_set)
    {
      requests.push_back(owner_cell_gid);
      requests.push_back(static_cast<uint64_t>(owner_face_index));
      requests.push_back(neighbor_cell_gid);
      requests.push_back(cg);
    }
  }

  const auto received_requests = MapAllToAll(pid_requests);

  std::map<int, std::vector<uint64_t>> pid_response_keys;
  std::map<int, std::vector<double>> pid_response_values;
  for (const auto& [pid, requests] : received_requests)
  {
    OpenSnLogicalErrorIf(requests.size() % 4 != 0, "Invalid CMFD face-current request buffer.");
    auto& response_keys = pid_response_keys[pid];
    auto& response_values = pid_response_values[pid];
    response_keys.reserve(requests.size());
    response_values.reserve(requests.size() / 4);

    for (std::size_t r = 0; r < requests.size(); r += 4)
    {
      const auto owner_cell_gid = requests[r];
      const auto requested_face_index = static_cast<std::size_t>(requests[r + 1]);
      const auto neighbor_cell_gid = requests[r + 2];
      const auto cg = static_cast<unsigned int>(requests[r + 3]);
      const auto& owner_cell = grid.cells[owner_cell_gid];
      OpenSnInvalidArgumentIf(owner_cell.partition_id != opensn::mpi_comm.rank(),
                              "CMFD face-current request references a nonlocal owner cell.");

      bool found_face = false;
      double outflow = 0.0;
      const auto& outflow_view = outflow_views[owner_cell.local_id];
      const auto wildcard_face_index = std::numeric_limits<std::size_t>::max();
      const std::size_t first_face =
        requested_face_index == wildcard_face_index ? 0 : requested_face_index;
      const std::size_t end_face = requested_face_index == wildcard_face_index
                                     ? owner_cell.faces.size()
                                     : requested_face_index + 1;
      for (std::size_t f = first_face; f < end_face; ++f)
      {
        const auto& face = owner_cell.faces[f];
        if ((face.has_neighbor and face.neighbor_id == neighbor_cell_gid) or
            (not face.has_neighbor and neighbor_cell_gid == std::numeric_limits<uint64_t>::max()))
        {
          for (unsigned int g = FineGroupBegin(cg); g < FineGroupEnd(cg); ++g)
            outflow += outflow_view.Get(f, first_group_ + g);
          found_face = true;
          break;
        }
      }

      OpenSnLogicalErrorIf(not found_face,
                           "Unable to satisfy CMFD face-current request for neighbor face.");
      response_keys.push_back(owner_cell_gid);
      response_keys.push_back(static_cast<uint64_t>(requested_face_index));
      response_keys.push_back(neighbor_cell_gid);
      response_keys.push_back(cg);
      response_values.push_back(outflow);
    }
  }

  const auto received_response_keys = MapAllToAll(pid_response_keys);
  const auto received_response_values = MapAllToAll(pid_response_values);
  face_current_cache_.clear();
  for (const auto& [pid, keys] : received_response_keys)
  {
    const auto values_it = received_response_values.find(pid);
    OpenSnLogicalErrorIf(values_it == received_response_values.end(),
                         "Missing CMFD face-current response values.");
    const auto& values = values_it->second;
    OpenSnLogicalErrorIf(keys.size() != 4 * values.size(),
                         "Invalid CMFD face-current response buffer.");

    for (std::size_t r = 0; r < values.size(); ++r)
    {
      const auto key_offset = 4 * r;
      const auto owner_cell_gid = keys[key_offset];
      const auto owner_face_index = static_cast<std::size_t>(keys[key_offset + 1]);
      const auto neighbor_cell_gid = keys[key_offset + 2];
      const auto cg = static_cast<unsigned int>(keys[key_offset + 3]);
      face_current_cache_[std::make_tuple(
        owner_cell_gid, owner_face_index, neighbor_cell_gid, cg)] = values[r];
    }
  }
}

void
CMFDAcceleration::AssembleOperator()
{
  const auto& xs_map = do_problem_.GetBlockID2XSMap();

  OpenSnPETScCall(MatZeroEntries(A_));

  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = *xs_map.at(coarse_cell.block_id);

    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      const auto row = MapDOF(coarse_cell.global_id, cg);
      double diag =
        CondensedRemovalXS(xs, coarse_phi_new_fine_, coarse_cell, cg) * coarse_cell.volume;

      for (std::size_t f = 0; f < coarse_cell.faces.size(); ++f)
      {
        const auto& face = coarse_cell.faces[f];
        if (face.has_neighbor)
        {
          const auto& neighbor_xs = *xs_map.at(face.neighbor_block_id);

          const double owner_distance_to_face =
            ProjectedDistanceToFace(coarse_cell.centroid, face.centroid, face.normal);
          const double neighbor_distance_to_face =
            ProjectedDistanceToFace(face.neighbor_centroid, face.centroid, face.normal);
          const std::vector<double> unweighted_spectrum;
          const double coeff = InterfaceDiffusionCoefficient(
            CondensedDiffusionCoefficient(xs, coarse_phi_new_fine_, coarse_cell, cg),
            CondensedDiffusionCoefficient(neighbor_xs, unweighted_spectrum, coarse_cell, cg),
            owner_distance_to_face,
            neighbor_distance_to_face,
            face.area);
          const double blend = std::clamp(current_closure_blend_, 0.0, 1.0);
          bool partial_current_closure_available = false;
          double partial_diag = 0.0;
          double partial_offdiag = 0.0;
          double current_correction = 0.0;
          if (not coarse_phi_new_.empty())
          {
            const auto neighbor_coarse_gid = face.neighbor_id;
            const auto phi_owner = coarse_phi_new_[coarse_cell.local_id * num_coarse_groups_ + cg];
            const auto phi_neighbor_it =
              coarse_phi_cache_.find(std::make_tuple(neighbor_coarse_gid, cg));
            if (phi_neighbor_it != coarse_phi_cache_.end())
            {
              const auto phi_neighbor = phi_neighbor_it->second;
              if (active_current_closure_ == "partial" and std::fabs(phi_owner) > 1.0e-14 and
                  std::fabs(phi_neighbor) > 1.0e-14)
              {
                const auto [owner_partial_current, neighbor_partial_current] =
                  ComputePartialOutwardCurrents(coarse_cell, f, cg);
                const double owner_coeff = owner_partial_current / phi_owner;
                const double neighbor_coeff = neighbor_partial_current / phi_neighbor;
                if (std::isfinite(owner_coeff) and std::isfinite(neighbor_coeff))
                {
                  partial_diag = owner_coeff;
                  partial_offdiag = -neighbor_coeff;
                  partial_current_closure_available = true;
                }
              }

              if ((active_current_closure_ == "net" or active_current_closure_ == "blend") and
                  blend < 1.0)
              {
                const auto phi_sum = phi_owner + phi_neighbor;
                if (std::fabs(phi_sum) > 1.0e-14)
                {
                  const double transport_current = ComputeOutwardCurrent(coarse_cell, f, cg);
                  current_correction =
                    (transport_current - coeff * (phi_owner - phi_neighbor)) / phi_sum;
                }
              }

              if (active_current_closure_ == "blend" and blend > 0.0 and
                  not partial_current_closure_available and std::fabs(phi_owner) > 1.0e-14 and
                  std::fabs(phi_neighbor) > 1.0e-14)
              {
                const auto [owner_partial_current, neighbor_partial_current] =
                  ComputePartialOutwardCurrents(coarse_cell, f, cg);
                const double owner_coeff = owner_partial_current / phi_owner;
                const double neighbor_coeff = neighbor_partial_current / phi_neighbor;
                if (std::isfinite(owner_coeff) and std::isfinite(neighbor_coeff))
                {
                  partial_diag = owner_coeff;
                  partial_offdiag = -neighbor_coeff;
                  partial_current_closure_available = true;
                }
              }
            }
          }

          const double net_diag = coeff + current_correction;
          const double net_offdiag = -coeff + current_correction;
          double face_diag = net_diag;
          double face_offdiag = net_offdiag;
          if (active_current_closure_ == "partial" and partial_current_closure_available)
          {
            face_diag = partial_diag;
            face_offdiag = partial_offdiag;
          }
          else if (active_current_closure_ == "blend" and partial_current_closure_available)
          {
            face_diag = (1.0 - blend) * net_diag + blend * partial_diag;
            face_offdiag = (1.0 - blend) * net_offdiag + blend * partial_offdiag;
          }

          diag += face_diag;
          const auto col = MapDOF(face.neighbor_id, cg);
          OpenSnPETScCall(MatSetValue(A_, row, col, face_offdiag, ADD_VALUES));
        }
        else
        {
          const auto bc_it = do_problem_.GetBoundaryDefinitions().find(face.neighbor_id);
          const bool reflective = bc_it != do_problem_.GetBoundaryDefinitions().end() and
                                  bc_it->second.type == LBSBoundaryType::REFLECTING;
          if (not reflective)
          {
            double boundary_coeff = 0.5 * face.area;
            if (not coarse_phi_new_.empty())
            {
              const auto phi_owner =
                coarse_phi_new_[coarse_cell.local_id * num_coarse_groups_ + cg];
              if (std::fabs(phi_owner) > 1.0e-14)
                boundary_coeff = ComputeOutwardCurrent(coarse_cell, f, cg) / phi_owner;
            }
            diag += boundary_coeff;
          }
        }
      }

      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        if (cgp != cg)
        {
          const auto sigma_s = CondensedScatterXS(xs, coarse_phi_new_fine_, coarse_cell, cg, cgp);
          const auto col = MapDOF(coarse_cell.global_id, cgp);
          OpenSnPETScCall(MatSetValue(A_, row, col, -sigma_s * coarse_cell.volume, ADD_VALUES));
        }

      OpenSnPETScCall(MatSetValue(A_, row, row, diag, ADD_VALUES));
    }
  }

  OpenSnPETScCall(MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY));
  OpenSnPETScCall(MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY));
  OpenSnPETScCall(KSPSetOperators(ksp_, A_, A_));
  OpenSnPETScCall(KSPSetUp(ksp_));
}

double
CMFDAcceleration::ComputeOutwardCurrent(const CMFDCoarseCell& coarse_cell,
                                        const std::size_t face_index,
                                        const unsigned int coarse_group) const
{
  const auto [owner_partial_current, neighbor_partial_current] =
    ComputePartialOutwardCurrents(coarse_cell, face_index, coarse_group);
  return owner_partial_current - neighbor_partial_current;
}

std::pair<double, double>
CMFDAcceleration::ComputePartialOutwardCurrents(const CMFDCoarseCell& coarse_cell,
                                                const std::size_t face_index,
                                                const unsigned int coarse_group) const
{
  const auto& grid = *do_problem_.GetGrid();
  const auto& coarse_face = coarse_cell.faces.at(face_index);
  const auto& outflow_views = do_problem_.GetCellOutflowViews();

  double owner_partial_current = 0.0;
  double neighbor_partial_current = 0.0;
  for (const auto& fine_face : coarse_face.fine_faces)
  {
    double owner_outflow = 0.0;
    if (fine_face.cell_partition_id == opensn::mpi_comm.rank())
    {
      const auto& fine_cell = grid.cells[fine_face.cell_id];
      for (unsigned int g = FineGroupBegin(coarse_group); g < FineGroupEnd(coarse_group); ++g)
        owner_outflow +=
          outflow_views[fine_cell.local_id].Get(fine_face.face_index, first_group_ + g);
    }
    else
    {
      const auto key =
        std::make_tuple(fine_face.cell_id,
                        fine_face.face_index,
                        fine_face.neighbor_id.value_or(std::numeric_limits<uint64_t>::max()),
                        coarse_group);
      const auto cached_outflow = face_current_cache_.find(key);
      OpenSnLogicalErrorIf(cached_outflow == face_current_cache_.end(),
                           "Unable to find matching remote owner face for CMFD current.");
      owner_outflow = cached_outflow->second;
    }
    owner_partial_current += owner_outflow;

    const auto neighbor_id = fine_face.neighbor_id;
    if (not neighbor_id.has_value())
    {
      continue;
    }

    if (fine_face.neighbor_partition_id == opensn::mpi_comm.rank())
    {
      const auto& neighbor_cell = grid.cells[*neighbor_id];
      bool found_face = false;
      for (std::size_t nf = 0; nf < neighbor_cell.faces.size(); ++nf)
      {
        const auto& neighbor_face = neighbor_cell.faces[nf];
        if (neighbor_face.has_neighbor and neighbor_face.neighbor_id == fine_face.cell_id)
        {
          double neighbor_outflow = 0.0;
          for (unsigned int g = FineGroupBegin(coarse_group); g < FineGroupEnd(coarse_group); ++g)
            neighbor_outflow += outflow_views[neighbor_cell.local_id].Get(nf, first_group_ + g);
          neighbor_partial_current += neighbor_outflow;
          found_face = true;
          break;
        }
      }
      OpenSnLogicalErrorIf(not found_face,
                           "Unable to find matching local neighbor face for CMFD current.");
    }
    else
    {
      const auto key = std::make_tuple(
        *neighbor_id, std::numeric_limits<std::size_t>::max(), fine_face.cell_id, coarse_group);
      const auto neighbor_outflow = face_current_cache_.find(key);
      OpenSnLogicalErrorIf(neighbor_outflow == face_current_cache_.end(),
                           "Unable to find matching remote neighbor face for CMFD current.");
      neighbor_partial_current += neighbor_outflow->second;
    }
  }

  return {owner_partial_current, neighbor_partial_current};
}

std::vector<double>
CMFDAcceleration::SolveCoarseSystem(const double k_eff,
                                    const std::vector<double>& coarse_source_phi) const
{
  OpenSnInvalidArgumentIf(coarse_source_phi.size() !=
                            coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse source scalar flux vector size mismatch.");

  AssembleRHS(rhs_, k_eff, coarse_source_phi);

  Vec x = nullptr;
  OpenSnPETScCall(VecDuplicate(rhs_, &x));
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      OpenSnPETScCall(VecSetValue(
        x, MapDOF(coarse_cell.global_id, cg), coarse_phi_new_[coarse_offset + cg], INSERT_VALUES));
  }
  OpenSnPETScCall(VecAssemblyBegin(x));
  OpenSnPETScCall(VecAssemblyEnd(x));

  KSPType ksp_type = nullptr;
  OpenSnPETScCall(KSPGetType(ksp_, &ksp_type));
  const bool use_nonzero_initial_guess =
    ksp_type == nullptr or std::string(ksp_type) != std::string(KSPPREONLY);
  OpenSnPETScCall(
    KSPSetInitialGuessNonzero(ksp_, use_nonzero_initial_guess ? PETSC_TRUE : PETSC_FALSE));
  OpenSnPETScCall(KSPSolve(ksp_, rhs_, x));
  PetscInt ksp_iterations = 0;
  OpenSnPETScCall(KSPGetIterationNumber(ksp_, &ksp_iterations));
  last_ksp_iterations_ = static_cast<int>(ksp_iterations);
  KSPConvergedReason ksp_reason = KSP_CONVERGED_ITERATING;
  OpenSnPETScCall(KSPGetConvergedReason(ksp_, &ksp_reason));
  last_ksp_converged_ = ksp_reason > 0;

  std::vector<int64_t> indices;
  indices.reserve(coarse_phi_new_.size());
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      indices.push_back(MapDOF(coarse_cell.global_id, cg));

  std::vector<double> coarse_phi;
  CopyGlobalVecToSTLvector(x, indices, coarse_phi);

  OpenSnPETScCall(VecDestroy(&x));
  return coarse_phi;
}

void
CMFDAcceleration::AssembleRHS(Vec rhs,
                              const double k_eff,
                              const std::vector<double>& coarse_source_phi) const
{
  const auto& xs_map = do_problem_.GetBlockID2XSMap();

  OpenSnInvalidArgumentIf(coarse_source_phi.size() !=
                            coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse source scalar flux vector size mismatch.");

  OpenSnPETScCall(VecSet(rhs, 0.0));
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = *xs_map.at(coarse_cell.block_id);
    if (not xs.IsFissionable())
      continue;

    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      double rhs_value = 0.0;
      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        rhs_value += CondensedProductionXS(xs, coarse_phi_new_fine_, coarse_cell, cg, cgp) *
                     coarse_source_phi[coarse_offset + cgp];
      rhs_value *= coarse_cell.volume / k_eff;
      const auto row = MapDOF(coarse_cell.global_id, cg);
      OpenSnPETScCall(VecSetValue(rhs, row, rhs_value, ADD_VALUES));
    }
  }
  OpenSnPETScCall(VecAssemblyBegin(rhs));
  OpenSnPETScCall(VecAssemblyEnd(rhs));
}

CMFDAcceleration::BalanceResidual
CMFDAcceleration::ComputeBalanceResidual(const double k_eff,
                                         const std::vector<double>& coarse_phi,
                                         const std::vector<double>& coarse_source_phi) const
{
  OpenSnInvalidArgumentIf(coarse_phi.size() != coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse scalar flux vector size mismatch.");

  Vec x = nullptr;
  Vec rhs = nullptr;
  Vec residual = nullptr;
  OpenSnPETScCall(VecDuplicate(rhs_, &x));
  OpenSnPETScCall(VecDuplicate(rhs_, &rhs));
  OpenSnPETScCall(VecDuplicate(rhs_, &residual));

  OpenSnPETScCall(VecSet(x, 0.0));
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      OpenSnPETScCall(VecSetValue(
        x, MapDOF(coarse_cell.global_id, cg), coarse_phi[coarse_offset + cg], INSERT_VALUES));
  }
  OpenSnPETScCall(VecAssemblyBegin(x));
  OpenSnPETScCall(VecAssemblyEnd(x));

  AssembleRHS(rhs, k_eff, coarse_source_phi);
  OpenSnPETScCall(MatMult(A_, x, residual));
  OpenSnPETScCall(VecAXPY(residual, -1.0, rhs));

  PetscReal max_abs = 0.0;
  PetscReal residual_l2 = 0.0;
  PetscReal rhs_l2 = 0.0;
  OpenSnPETScCall(VecNorm(residual, NORM_INFINITY, &max_abs));
  OpenSnPETScCall(VecNorm(residual, NORM_2, &residual_l2));
  OpenSnPETScCall(VecNorm(rhs, NORM_2, &rhs_l2));

  OpenSnPETScCall(VecDestroy(&residual));
  OpenSnPETScCall(VecDestroy(&rhs));
  OpenSnPETScCall(VecDestroy(&x));

  return {static_cast<double>(max_abs),
          static_cast<double>(residual_l2 / std::max<PetscReal>(rhs_l2, 1.0e-30))};
}

CMFDAcceleration::BalanceResidual
CMFDAcceleration::ComputeTransportCurrentBalanceResidual(
  const double k_eff,
  const std::vector<double>& coarse_phi,
  const std::vector<double>& coarse_source_phi) const
{
  OpenSnInvalidArgumentIf(coarse_phi.size() != coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse scalar flux vector size mismatch.");
  OpenSnInvalidArgumentIf(coarse_source_phi.size() !=
                            coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse source scalar flux vector size mismatch.");

  const auto& xs_map = do_problem_.GetBlockID2XSMap();
  double local_max_abs = 0.0;
  double local_residual_l2 = 0.0;
  double local_rhs_l2 = 0.0;

  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = *xs_map.at(coarse_cell.block_id);
    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;

    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      const double phi = coarse_phi[coarse_offset + cg];
      double lhs =
        CondensedRemovalXS(xs, coarse_phi_new_fine_, coarse_cell, cg) * coarse_cell.volume * phi;

      for (std::size_t f = 0; f < coarse_cell.faces.size(); ++f)
      {
        const auto& face = coarse_cell.faces[f];
        if (face.has_neighbor)
          lhs += ComputeOutwardCurrent(coarse_cell, f, cg);
        else
        {
          const auto bc_it = do_problem_.GetBoundaryDefinitions().find(face.neighbor_id);
          const bool reflective = bc_it != do_problem_.GetBoundaryDefinitions().end() and
                                  bc_it->second.type == LBSBoundaryType::REFLECTING;
          if (not reflective)
            lhs += ComputeOutwardCurrent(coarse_cell, f, cg);
        }
      }

      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        if (cgp != cg)
        {
          const auto sigma_s = CondensedScatterXS(xs, coarse_phi_new_fine_, coarse_cell, cg, cgp);
          lhs -= sigma_s * coarse_cell.volume * coarse_phi[coarse_offset + cgp];
        }

      double rhs = 0.0;
      if (xs.IsFissionable())
        for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
          rhs += CondensedProductionXS(xs, coarse_phi_new_fine_, coarse_cell, cg, cgp) *
                 coarse_source_phi[coarse_offset + cgp];
      rhs *= coarse_cell.volume / k_eff;

      const double residual = lhs - rhs;
      local_max_abs = std::max(local_max_abs, std::fabs(residual));
      local_residual_l2 += residual * residual;
      local_rhs_l2 += rhs * rhs;
    }
  }

  double global_max_abs = 0.0;
  double global_residual_l2 = 0.0;
  double global_rhs_l2 = 0.0;
  mpi_comm.all_reduce(local_max_abs, global_max_abs, mpi::op::max<double>());
  mpi_comm.all_reduce(local_residual_l2, global_residual_l2, mpi::op::sum<double>());
  mpi_comm.all_reduce(local_rhs_l2, global_rhs_l2, mpi::op::sum<double>());

  return {global_max_abs, std::sqrt(global_residual_l2 / std::max(global_rhs_l2, 1.0e-60))};
}

double
CMFDAcceleration::ComputeCoarseFissionProduction(const std::vector<double>& coarse_phi) const
{
  const auto& xs_map = do_problem_.GetBlockID2XSMap();

  OpenSnInvalidArgumentIf(coarse_phi.size() != coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse scalar flux vector size mismatch.");

  double local_production = 0.0;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = *xs_map.at(coarse_cell.block_id);
    if (not xs.IsFissionable())
      continue;

    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        local_production += CondensedProductionXS(xs, coarse_phi_new_fine_, coarse_cell, cg, cgp) *
                            coarse_phi[coarse_offset + cgp] * coarse_cell.volume;
  }

  double global_production = 0.0;
  mpi_comm.all_reduce(local_production, global_production, mpi::op::sum<double>());
  return global_production;
}

CMFDAcceleration::FluxUpdateDiagnostics
CMFDAcceleration::AnalyzeFluxUpdate(const std::vector<double>& phi, const double k_eff) const
{
  const auto& transport_views = do_problem_.GetCellTransportViews();
  const auto& grid = *do_problem_.GetGrid();
  double local_min_phi = std::numeric_limits<double>::max();
  bool local_nonfinite = false;
  for (const auto& cell : grid.local_cells)
  {
    const auto& transport_view = transport_views[cell.local_id];
    for (int i = 0; i < transport_view.GetNumNodes(); ++i)
    {
      const auto phi_map = transport_view.MapDOF(i, 0, first_group_);
      for (unsigned int g = 0; g < num_groups_; ++g)
      {
        const auto value = phi[phi_map + g];
        local_nonfinite = local_nonfinite or not std::isfinite(value);
        local_min_phi = std::min(local_min_phi, value);
      }
    }
  }

  FluxUpdateDiagnostics diagnostics;
  int global_nonfinite = 0;
  int global_invalid_k = 0;
  opensn::mpi_comm.all_reduce(local_min_phi, diagnostics.min_scalar_flux, mpi::op::min<double>());
  opensn::mpi_comm.all_reduce(local_nonfinite ? 1 : 0, global_nonfinite, mpi::op::max<int>());
  opensn::mpi_comm.all_reduce(
    k_eff <= 0.0 or not std::isfinite(k_eff) ? 1 : 0, global_invalid_k, mpi::op::max<int>());
  diagnostics.has_nonfinite_flux = global_nonfinite != 0;
  diagnostics.has_invalid_k = global_invalid_k != 0;
  return diagnostics;
}

double
CMFDAcceleration::ApplyFluxCorrectionWithDamping(const std::vector<double>& transport_phi_new,
                                                 const std::vector<double>& coarse_phi,
                                                 const double old_fission_production,
                                                 const double transport_k_eff,
                                                 double& applied_damping,
                                                 unsigned int& attempts,
                                                 FluxUpdateDiagnostics& diagnostics,
                                                 bool& skipped_correction,
                                                 std::string& skip_reason)
{
  std::vector<double> unrelaxed_coarse_delta_phi(coarse_phi.size(), 0.0);
  for (size_t i = 0; i < coarse_phi.size(); ++i)
    unrelaxed_coarse_delta_phi[i] = coarse_phi[i] - coarse_phi_new_[i];

  double damping = relaxation_;
  attempts = 0;
  skipped_correction = false;
  skip_reason.clear();
  applied_damping = 0.0;
  std::string rejection_reason = "correction_damping_exhausted";
  const auto transport_diagnostics = AnalyzeFluxUpdate(transport_phi_new, transport_k_eff);
  OpenSnLogicalErrorIf(transport_diagnostics.has_invalid_k,
                       "CMFD input transport update has an invalid k-eigenvalue.");
  OpenSnLogicalErrorIf(transport_diagnostics.has_nonfinite_flux,
                       "CMFD input transport update has non-finite scalar flux values.");
  const double min_allowed_scalar_flux = std::min(
    -negative_flux_tolerance_, transport_diagnostics.min_scalar_flux - negative_flux_tolerance_);

  for (; attempts < correction_max_attempts_; ++attempts)
  {
    if (attempts > 0 and damping < correction_min_damping_)
      break;

    std::vector<double> coarse_delta_phi(unrelaxed_coarse_delta_phi.size(), 0.0);
    for (size_t i = 0; i < unrelaxed_coarse_delta_phi.size(); ++i)
      coarse_delta_phi[i] = damping * unrelaxed_coarse_delta_phi[i];
    std::vector<double> candidate_coarse_phi(coarse_phi_new_.size(), 0.0);
    for (size_t i = 0; i < candidate_coarse_phi.size(); ++i)
      candidate_coarse_phi[i] = coarse_phi_new_[i] + coarse_delta_phi[i];

    auto candidate_phi_new = transport_phi_new;
    CMFDProlongateScalarFluxRatio(do_problem_,
                                  first_group_,
                                  num_groups_,
                                  group_aggregation_size_,
                                  coarse_mesh_,
                                  coarse_phi_new_,
                                  candidate_coarse_phi,
                                  candidate_phi_new);
    const double candidate_fission_production =
      ComputeFissionProduction(do_problem_, candidate_phi_new);
    const double candidate_k_eff =
      candidate_fission_production / old_fission_production * solver_->GetEigenvalue();
    diagnostics = AnalyzeFluxUpdate(candidate_phi_new, candidate_k_eff);
    if (not diagnostics.has_invalid_k and not diagnostics.has_nonfinite_flux and
        diagnostics.min_scalar_flux >= min_allowed_scalar_flux)
    {
      phi_new_local_ = std::move(candidate_phi_new);
      applied_damping = damping;
      ++attempts;
      return candidate_k_eff;
    }

    if (diagnostics.has_invalid_k)
      rejection_reason = "invalid_k";
    else if (diagnostics.has_nonfinite_flux)
      rejection_reason = "nonfinite_flux";
    else
      rejection_reason = "negative_flux_guard";

    damping *= 0.5;
  }

  diagnostics = transport_diagnostics;
  phi_new_local_ = transport_phi_new;
  skipped_correction = true;
  skip_reason = rejection_reason;
  if (verbose_ or outer_iteration_ < 5 or outer_iteration_ % 25 == 0)
    log.Log0Warning() << no_wrap << "CMFD correction rejected after " << attempts
                      << " damping attempts; using the unaccelerated transport update. "
                      << "Transport minimum scalar flux = " << diagnostics.min_scalar_flux << ".";
  return transport_k_eff;
}

} // namespace opensn
