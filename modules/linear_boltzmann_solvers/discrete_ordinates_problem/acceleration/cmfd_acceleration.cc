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
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
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

} // namespace

OpenSnRegisterObjectInNamespace(lbs, CMFDAcceleration);

InputParameters
CMFDAcceleration::GetInputParameters()
{
  auto params = DiscreteOrdinatesKEigenAcceleration::GetInputParameters();

  params.ChangeExistingParamToOptional("pi_max_its", 5);

  params.AddOptionalParameter(
    "coarse_mesh", "local_aggregation", "CMFD coarse-mesh construction method.");
  params.AddOptionalParameter(
    "aggregation_size", 16, "Target number of fine cells per aggregated CMFD coarse cell.");
  params.AddOptionalParameter(
    "relaxation", 1.0, "Relaxation factor applied to the CMFD scalar-flux correction.");
  params.AddOptionalParameter("correction_max_attempts",
                              10,
                              "Maximum CMFD correction damping attempts before skipping the "
                              "correction for the current transport update.");
  params.AddOptionalParameter(
    "correction_min_damping",
    1.0e-4,
    "Minimum CMFD correction damping factor considered during correction limiting.");
  params.AddOptionalParameter("negative_flux_tolerance",
                              1.0e-6,
                              "Allowed scalar-flux undershoot for accepting a CMFD correction.");
  params.AddOptionalParameter(
    "adaptive_relaxation", true, "If true, adapt the starting CMFD correction relaxation.");
  params.AddOptionalParameter(
    "adaptive_relaxation_min", 0.25, "Minimum starting relaxation used by adaptive relaxation.");
  params.AddOptionalParameter(
    "adaptive_relaxation_max", 1.0, "Maximum starting relaxation used by adaptive relaxation.");
  params.AddOptionalParameter(
    "adaptive_relaxation_growth",
    1.25,
    "Multiplier used to increase adaptive relaxation after repeated strong accepted corrections.");
  params.AddOptionalParameter("adaptive_relaxation_reduction",
                              0.75,
                              "Multiplier used to decrease adaptive relaxation after rejected or "
                              "strongly damped corrections.");
  params.AddOptionalParameter(
    "adaptive_relaxation_accept_fraction",
    0.5,
    "Fraction of the starting relaxation required to treat an accepted correction as strong.");
  params.AddOptionalParameter(
    "adaptive_relaxation_successes_to_grow",
    2,
    "Number of consecutive strong accepted corrections before increasing adaptive relaxation.");
  params.AddOptionalParameter("inactive_iterations",
                              1,
                              "Number of initial power iterations before applying CMFD "
                              "corrections. Transport updates are still used.");
  params.AddOptionalParameter("update_scheme",
                              true,
                              "If true, configures the transport solve as a loose update solve "
                              "and applies CMFD after each update.");
  params.AddOptionalParameter(
    "update_wgs_max_its", 4, "Maximum WGS iterations used when update_scheme is true.");
  params.AddOptionalParameter(
    "update_wgs_abs_tol", 1.0e-4, "WGS absolute tolerance used when update_scheme is true.");
  params.AddOptionalParameter(
    "coarse_solver_policy",
    "auto",
    "CMFD coarse solver policy: auto, direct, iterative, or petsc_options.");
  params.AddOptionalParameter("direct_coarse_solve_threshold",
                              5000,
                              "Maximum global CMFD unknown count for auto direct coarse solves.");
  params.ConstrainParameterRange("coarse_mesh",
                                 AllowableRangeList::New({"identity", "local_aggregation"}));
  params.ConstrainParameterRange(
    "coarse_solver_policy",
    AllowableRangeList::New({"auto", "direct", "iterative", "petsc_options"}));
  params.ConstrainParameterRange("aggregation_size", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("relaxation", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("correction_max_attempts", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("correction_min_damping", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("negative_flux_tolerance", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("adaptive_relaxation_min", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("adaptive_relaxation_max", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("adaptive_relaxation_growth", AllowableRangeLowLimit::New(1.0));
  params.ConstrainParameterRange("adaptive_relaxation_reduction", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("adaptive_relaxation_accept_fraction",
                                 AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("adaptive_relaxation_successes_to_grow",
                                 AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("inactive_iterations", AllowableRangeLowLimit::New(0));
  params.ConstrainParameterRange("update_wgs_max_its", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("update_wgs_abs_tol", AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("direct_coarse_solve_threshold", AllowableRangeLowLimit::New(1));
  params.ChangeExistingParamToOptional("name", "CMFDAcceleration");

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
    aggregation_size_(params.GetParamValue<int>("aggregation_size")),
    relaxation_(params.GetParamValue<double>("relaxation")),
    correction_max_attempts_(params.GetParamValue<unsigned int>("correction_max_attempts")),
    correction_min_damping_(params.GetParamValue<double>("correction_min_damping")),
    negative_flux_tolerance_(params.GetParamValue<double>("negative_flux_tolerance")),
    adaptive_relaxation_(params.GetParamValue<bool>("adaptive_relaxation")),
    adaptive_relaxation_min_(params.GetParamValue<double>("adaptive_relaxation_min")),
    adaptive_relaxation_max_(params.GetParamValue<double>("adaptive_relaxation_max")),
    adaptive_relaxation_growth_(params.GetParamValue<double>("adaptive_relaxation_growth")),
    adaptive_relaxation_reduction_(params.GetParamValue<double>("adaptive_relaxation_reduction")),
    adaptive_relaxation_accept_fraction_(
      params.GetParamValue<double>("adaptive_relaxation_accept_fraction")),
    adaptive_relaxation_successes_to_grow_(
      params.GetParamValue<unsigned int>("adaptive_relaxation_successes_to_grow")),
    inactive_iterations_(params.GetParamValue<unsigned int>("inactive_iterations")),
    update_scheme_(params.GetParamValue<bool>("update_scheme")),
    update_wgs_max_its_(params.GetParamValue<unsigned int>("update_wgs_max_its")),
    update_wgs_abs_tol_(params.GetParamValue<double>("update_wgs_abs_tol")),
    coarse_solver_policy_(params.GetParamValue<std::string>("coarse_solver_policy")),
    direct_coarse_solve_threshold_(
      static_cast<std::size_t>(params.GetParamValue<int>("direct_coarse_solve_threshold"))),
    num_groups_(do_problem_.GetNumGroups())
{
  OpenSnInvalidArgumentIf(adaptive_relaxation_max_ < adaptive_relaxation_min_,
                          "CMFD adaptive_relaxation_max must be >= adaptive_relaxation_min.");
  OpenSnInvalidArgumentIf(adaptive_relaxation_reduction_ > 1.0,
                          "CMFD adaptive_relaxation_reduction must be <= 1.0.");
  OpenSnInvalidArgumentIf(adaptive_relaxation_accept_fraction_ > 1.0,
                          "CMFD adaptive_relaxation_accept_fraction must be <= 1.0.");
  current_relaxation_ = std::clamp(relaxation_, adaptive_relaxation_min_, adaptive_relaxation_max_);
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
  if (coarse_mesh_type_ == "identity")
    coarse_mesh_ = CMFDCoarseMesh::BuildIdentity(*do_problem_.GetGrid());
  else if (coarse_mesh_type_ == "local_aggregation")
    coarse_mesh_ = CMFDCoarseMesh::BuildLocalAggregation(*do_problem_.GetGrid(), aggregation_size_);

  ApplyTransportUpdateScheme();
  InitializeLinearSystem();
  AssembleOperator();

  log.Log() << "Initialized CMFD coarse mesh with " << coarse_mesh_.NumLocalCells()
            << " local coarse cells.";
  LogCoarseMeshDiagnostics(ComputeCoarseMeshDiagnostics());
}

void
CMFDAcceleration::PreExecute()
{
  coarse_phi_old_.clear();
  coarse_phi_new_.clear();
  coarse_phi_cache_.clear();
  face_current_cache_.clear();
  outer_iteration_ = 0;
  last_restrict_old_time_ = 0.0;
  transport_start_time_ = 0.0;
}

void
CMFDAcceleration::PrePowerIteration()
{
  const double t0 = NowSeconds();
  coarse_phi_old_ =
    CMFDRestrictScalarFlux(do_problem_, first_group_, num_groups_, coarse_mesh_, phi_old_local_);
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
  coarse_phi_new_ =
    CMFDRestrictScalarFlux(do_problem_, first_group_, num_groups_, coarse_mesh_, phi_new_local_);
  const double restrict_new_time = NowSeconds() - t0;

  t0 = NowSeconds();
  BuildCoarseFluxCache();
  const double coarse_flux_cache_time = NowSeconds() - t0;

  t0 = NowSeconds();
  BuildFaceCurrentCache();
  const double face_current_cache_time = NowSeconds() - t0;

  t0 = NowSeconds();
  AssembleOperator();
  const double assemble_time = NowSeconds() - t0;

  t0 = NowSeconds();
  const auto transport_balance_residual =
    ComputeBalanceResidual(solver_->GetEigenvalue(), coarse_phi_new_, coarse_phi_old_);
  const double residual_time = NowSeconds() - t0;
  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
    log.Log() << "CMFD restricted transport balance residual: max_abs = "
              << transport_balance_residual.max_abs
              << ", relative_l2 = " << transport_balance_residual.relative_l2;

  t0 = NowSeconds();
  const double F_old = ComputeFissionProduction(do_problem_, phi_old_local_);
  const double F_transport_new = ComputeFissionProduction(do_problem_, phi_new_local_);
  OpenSnLogicalErrorIf(F_old == 0.0, "CMFDAcceleration cannot update k with zero old production.");
  const double transport_k_eff = F_transport_new / F_old * solver_->GetEigenvalue();
  const double fission_time = NowSeconds() - t0;

  if (outer_iteration_ < inactive_iterations_)
  {
    for (const auto& groupset : groupsets_)
      LBSVecOps::GSScopedCopyPrimarySTLvectors(
        do_problem_, groupset, phi_new_local_, phi_old_local_);
    last_update_allows_convergence_ = false;

    if (verbose_ and outer_iteration_ < 5)
      log.Log() << "CMFD inactive during warmup iteration " << outer_iteration_ << " of "
                << inactive_iterations_ << ".";

    ++outer_iteration_;
    return transport_k_eff;
  }

  auto coarse_phi_m = coarse_phi_new_;
  auto coarse_phi = coarse_phi_new_;
  double k_eff = solver_->GetEigenvalue();
  double coarse_production_old = ComputeCoarseFissionProduction(coarse_phi_m);
  OpenSnLogicalErrorIf(coarse_production_old == 0.0,
                       "CMFDAcceleration cannot iterate with zero coarse production.");
  int coarse_pi_iterations = 0;
  int coarse_ksp_iterations = 0;
  t0 = NowSeconds();
  for (int m = 0; m < pi_max_its_; ++m)
  {
    coarse_phi = SolveCoarseSystem(k_eff, coarse_phi_m);
    coarse_ksp_iterations += last_ksp_iterations_;
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

  t0 = NowSeconds();
  const double starting_relaxation = adaptive_relaxation_ ? current_relaxation_ : relaxation_;
  double applied_damping = 0.0;
  unsigned int correction_attempts = 0;
  FluxUpdateDiagnostics correction_diagnostics;
  bool skipped_correction = false;
  k_eff = ApplyFluxCorrectionWithDamping(phi_new_local_,
                                         coarse_phi,
                                         transport_k_eff,
                                         k_eff,
                                         applied_damping,
                                         correction_attempts,
                                         correction_diagnostics,
                                         skipped_correction);
  last_update_allows_convergence_ = not skipped_correction;
  UpdateAdaptiveRelaxation(starting_relaxation, applied_damping, skipped_correction);
  for (const auto& groupset : groupsets_)
    LBSVecOps::GSScopedCopyPrimarySTLvectors(do_problem_, groupset, phi_new_local_, phi_old_local_);
  const double update_time = NowSeconds() - t0;
  const double post_time = NowSeconds() - post_start_time;

  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
    LogCorrectionDiagnostics(starting_relaxation,
                             applied_damping,
                             correction_attempts,
                             correction_diagnostics,
                             skipped_correction);

  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
    LogTimingSummary(transport_time,
                     last_restrict_old_time_,
                     restrict_new_time,
                     coarse_flux_cache_time,
                     face_current_cache_time,
                     assemble_time,
                     residual_time,
                     fission_time,
                     coarse_pi_time,
                     update_time,
                     post_time,
                     coarse_pi_iterations,
                     coarse_ksp_iterations,
                     transport_balance_residual);

  ++outer_iteration_;
  return k_eff;
}

std::size_t
CMFDAcceleration::LocalUnknownCount() const
{
  return coarse_mesh_.NumLocalCells() * num_groups_;
}

std::size_t
CMFDAcceleration::GlobalUnknownCount() const
{
  return coarse_mesh_.NumGlobalCells() * num_groups_;
}

PetscInt
CMFDAcceleration::MapDOF(const uint64_t coarse_cell_global_id, const unsigned int group) const
{
  return static_cast<PetscInt>(coarse_cell_global_id * num_groups_ + group);
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
  PC pc = nullptr;
  OpenSnPETScCall(KSPGetPC(ksp_, &pc));
  const bool direct_coarse_solve =
    coarse_solver_policy_ == "direct" or
    (coarse_solver_policy_ == "auto" and GlobalUnknownCount() <= direct_coarse_solve_threshold_);
  const bool petsc_options_policy = coarse_solver_policy_ == "petsc_options";
  if (direct_coarse_solve)
  {
    OpenSnPETScCall(KSPSetType(ksp_, KSPPREONLY));
    OpenSnPETScCall(PCSetType(pc, PCLU));
  }
  else
  {
    OpenSnPETScCall(KSPSetType(ksp_, KSPGMRES));
    OpenSnPETScCall(PCSetType(pc, PCJACOBI));
  }
  OpenSnPETScCall(KSPSetTolerances(ksp_, l_abs_tol_, PETSC_DEFAULT, PETSC_DEFAULT, max_iters_));
  if ((petsc_options_policy or petsc_options_ != "ssss") and not petsc_options_.empty())
    OpenSnPETScCall(PetscOptionsInsertString(nullptr, petsc_options_.c_str()));
  OpenSnPETScCall(KSPSetFromOptions(ksp_));
}

void
CMFDAcceleration::ApplyTransportUpdateScheme()
{
  if (not update_scheme_)
    return;

  for (std::size_t gsid = 0; gsid < do_problem_.GetNumGroupsets(); ++gsid)
  {
    auto& groupset = do_problem_.GetGroupset(gsid);
    groupset.max_iterations = update_wgs_max_its_;
    groupset.residual_tolerance = update_wgs_abs_tol_;

    auto wgs_solver =
      std::dynamic_pointer_cast<PETScLinearSolver>(do_problem_.GetWGSSolver(groupset.id));
    OpenSnLogicalErrorIf(not wgs_solver, "CMFD update scheme requires a PETSc WGS solver.");
    auto& tolerance_options = wgs_solver->GetToleranceOptions();
    tolerance_options.maximum_iterations = static_cast<PetscInt>(update_wgs_max_its_);
    tolerance_options.residual_absolute = update_wgs_abs_tol_;
  }

  log.Log() << "CMFD update scheme configured " << do_problem_.GetNumGroupsets()
            << " groupset WGS solvers with max iterations = " << update_wgs_max_its_
            << ", absolute tolerance = " << update_wgs_abs_tol_ << ".";
}

void
CMFDAcceleration::BuildCoarseFluxCache()
{
  coarse_phi_cache_.clear();

  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto offset = coarse_cell.local_id * num_groups_;
    for (unsigned int g = 0; g < num_groups_; ++g)
      coarse_phi_cache_[std::make_tuple(coarse_cell.global_id, g)] = coarse_phi_new_[offset + g];
  }

  std::map<int, std::set<std::pair<uint64_t, unsigned int>>> pid_request_sets;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    for (const auto& face : coarse_cell.faces)
    {
      if (not face.has_neighbor or face.neighbor_partition_id == opensn::mpi_comm.rank())
        continue;

      auto& requests = pid_request_sets[face.neighbor_partition_id];
      for (unsigned int g = 0; g < num_groups_; ++g)
        requests.emplace(face.neighbor_id, g);
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
      const auto g = static_cast<unsigned int>(requests[r + 1]);
      OpenSnInvalidArgumentIf(not coarse_mesh_.HasLocalCoarseCell(coarse_cell_gid),
                              "CMFD coarse-flux request references a nonlocal coarse cell.");

      const auto coarse_local_id = coarse_mesh_.LocalCellFromGlobalID(coarse_cell_gid).local_id;
      const auto value = coarse_phi_new_[coarse_local_id * num_groups_ + g];
      response_keys.push_back(coarse_cell_gid);
      response_keys.push_back(g);
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
      const auto g = static_cast<unsigned int>(keys[key_offset + 1]);
      coarse_phi_cache_[std::make_tuple(coarse_cell_gid, g)] = values[r];
    }
  }
}

void
CMFDAcceleration::BuildFaceCurrentCache()
{
  const auto& grid = *do_problem_.GetGrid();
  const auto& outflow_views = do_problem_.GetCellOutflowViews();

  std::map<int, std::set<std::tuple<uint64_t, uint64_t, unsigned int>>> pid_request_sets;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    for (const auto& face : coarse_cell.faces)
    {
      if (not face.has_neighbor)
        continue;

      for (const auto& fine_face : face.fine_faces)
      {
        if (not fine_face.neighbor_id.has_value())
          continue;
        const uint64_t neighbor_id = fine_face.neighbor_id.value();
        if (grid.IsCellLocal(neighbor_id))
          continue;

        const auto& neighbor_cell = grid.cells[neighbor_id];
        auto& requests = pid_request_sets[neighbor_cell.partition_id];
        for (unsigned int g = 0; g < num_groups_; ++g)
          requests.emplace(neighbor_id, fine_face.cell_id, g);
      }
    }
  }

  std::map<int, std::vector<uint64_t>> pid_requests;
  for (const auto& [pid, request_set] : pid_request_sets)
  {
    auto& requests = pid_requests[pid];
    requests.reserve(3 * request_set.size());
    for (const auto& [owner_cell_gid, neighbor_cell_gid, g] : request_set)
    {
      requests.push_back(owner_cell_gid);
      requests.push_back(neighbor_cell_gid);
      requests.push_back(g);
    }
  }

  const auto received_requests = MapAllToAll(pid_requests);

  std::map<int, std::vector<uint64_t>> pid_response_keys;
  std::map<int, std::vector<double>> pid_response_values;
  for (const auto& [pid, requests] : received_requests)
  {
    OpenSnLogicalErrorIf(requests.size() % 3 != 0, "Invalid CMFD face-current request buffer.");
    auto& response_keys = pid_response_keys[pid];
    auto& response_values = pid_response_values[pid];
    response_keys.reserve(requests.size());
    response_values.reserve(requests.size() / 3);

    for (std::size_t r = 0; r < requests.size(); r += 3)
    {
      const auto owner_cell_gid = requests[r];
      const auto neighbor_cell_gid = requests[r + 1];
      const auto g = static_cast<unsigned int>(requests[r + 2]);
      const auto& owner_cell = grid.cells[owner_cell_gid];
      OpenSnInvalidArgumentIf(owner_cell.partition_id != opensn::mpi_comm.rank(),
                              "CMFD face-current request references a nonlocal owner cell.");

      bool found_face = false;
      double outflow = 0.0;
      const auto& outflow_view = outflow_views[owner_cell.local_id];
      for (std::size_t f = 0; f < owner_cell.faces.size(); ++f)
      {
        const auto& face = owner_cell.faces[f];
        if (face.has_neighbor and face.neighbor_id == neighbor_cell_gid)
        {
          outflow = outflow_view.Get(f, g);
          found_face = true;
          break;
        }
      }

      OpenSnLogicalErrorIf(not found_face,
                           "Unable to satisfy CMFD face-current request for neighbor face.");
      response_keys.push_back(owner_cell_gid);
      response_keys.push_back(neighbor_cell_gid);
      response_keys.push_back(g);
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
    OpenSnLogicalErrorIf(keys.size() != 3 * values.size(),
                         "Invalid CMFD face-current response buffer.");

    for (std::size_t r = 0; r < values.size(); ++r)
    {
      const auto key_offset = 3 * r;
      const auto owner_cell_gid = keys[key_offset];
      const auto neighbor_cell_gid = keys[key_offset + 1];
      const auto g = static_cast<unsigned int>(keys[key_offset + 2]);
      face_current_cache_[std::make_tuple(owner_cell_gid, neighbor_cell_gid, g)] = values[r];
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
    const auto& diffusion_coeff = xs.GetDiffusionCoefficient();
    const auto& sigma_removal = xs.GetSigmaRemoval();
    const auto& scatter = xs.GetTransferMatrix(0);

    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      const auto row = MapDOF(coarse_cell.global_id, g);
      double diag = sigma_removal[g] * coarse_cell.volume;

      for (std::size_t f = 0; f < coarse_cell.faces.size(); ++f)
      {
        const auto& face = coarse_cell.faces[f];
        if (face.has_neighbor)
        {
          const auto& neighbor_xs = *xs_map.at(face.neighbor_block_id);
          const auto& neighbor_diffusion_coeff = neighbor_xs.GetDiffusionCoefficient();

          const double owner_distance_to_face =
            ProjectedDistanceToFace(coarse_cell.centroid, face.centroid, face.normal);
          const double neighbor_distance_to_face =
            ProjectedDistanceToFace(face.neighbor_centroid, face.centroid, face.normal);
          const double coeff = InterfaceDiffusionCoefficient(diffusion_coeff[g],
                                                             neighbor_diffusion_coeff[g],
                                                             owner_distance_to_face,
                                                             neighbor_distance_to_face,
                                                             face.area);
          double current_correction = 0.0;
          if (not coarse_phi_new_.empty())
          {
            const auto neighbor_coarse_gid = face.neighbor_id;
            const auto phi_owner = coarse_phi_new_[coarse_cell.local_id * num_groups_ + g];
            const auto phi_neighbor_it =
              coarse_phi_cache_.find(std::make_tuple(neighbor_coarse_gid, g));
            if (phi_neighbor_it != coarse_phi_cache_.end())
            {
              const auto phi_neighbor = phi_neighbor_it->second;
              const auto phi_sum = phi_owner + phi_neighbor;
              if (std::fabs(phi_sum) > 1.0e-14)
              {
                const double transport_current = ComputeOutwardCurrent(coarse_cell, f, g);
                current_correction =
                  (transport_current - coeff * (phi_owner - phi_neighbor)) / phi_sum;
              }
            }
          }

          diag += coeff + current_correction;
          const auto col = MapDOF(face.neighbor_id, g);
          OpenSnPETScCall(MatSetValue(A_, row, col, -coeff + current_correction, ADD_VALUES));
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
              const auto phi_owner = coarse_phi_new_[coarse_cell.local_id * num_groups_ + g];
              if (std::fabs(phi_owner) > 1.0e-14)
                boundary_coeff = ComputeOutwardCurrent(coarse_cell, f, g) / phi_owner;
            }
            diag += boundary_coeff;
          }
        }
      }

      for (const auto& [_, gp, sigma_s] : scatter.Row(g))
        if (gp >= first_group_ and gp < first_group_ + num_groups_ and gp != g)
        {
          const auto col = MapDOF(coarse_cell.global_id, gp - first_group_);
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
                                        const unsigned int group) const
{
  const auto& grid = *do_problem_.GetGrid();
  const auto& coarse_face = coarse_cell.faces.at(face_index);
  const auto& outflow_views = do_problem_.GetCellOutflowViews();
  const auto g = first_group_ + group;

  double current = 0.0;
  for (const auto& fine_face : coarse_face.fine_faces)
  {
    const auto& fine_cell = grid.cells[fine_face.cell_id];
    const double owner_outflow = outflow_views[fine_cell.local_id].Get(fine_face.face_index, g);

    if (not fine_face.neighbor_id.has_value())
    {
      current += owner_outflow;
      continue;
    }

    const uint64_t neighbor_id = fine_face.neighbor_id.value();
    if (grid.IsCellLocal(neighbor_id))
    {
      const auto& neighbor_cell = grid.cells[neighbor_id];
      bool found_face = false;
      for (std::size_t nf = 0; nf < neighbor_cell.faces.size(); ++nf)
      {
        const auto& neighbor_face = neighbor_cell.faces[nf];
        if (neighbor_face.has_neighbor and neighbor_face.neighbor_id == fine_cell.global_id)
        {
          current += owner_outflow - outflow_views[neighbor_cell.local_id].Get(nf, g);
          found_face = true;
          break;
        }
      }
      OpenSnLogicalErrorIf(not found_face,
                           "Unable to find matching local neighbor face for CMFD current.");
    }
    else
    {
      const auto key = std::make_tuple(neighbor_id, fine_cell.global_id, g);
      const auto neighbor_outflow = face_current_cache_.find(key);
      OpenSnLogicalErrorIf(neighbor_outflow == face_current_cache_.end(),
                           "Unable to find matching remote neighbor face for CMFD current.");
      current += owner_outflow - neighbor_outflow->second;
    }
  }

  return current;
}

std::vector<double>
CMFDAcceleration::SolveCoarseSystem(const double k_eff,
                                    const std::vector<double>& coarse_source_phi) const
{
  OpenSnInvalidArgumentIf(coarse_source_phi.size() != coarse_mesh_.NumLocalCells() * num_groups_,
                          "Coarse source scalar flux vector size mismatch.");

  AssembleRHS(rhs_, k_eff, coarse_source_phi);

  Vec x = nullptr;
  OpenSnPETScCall(VecDuplicate(rhs_, &x));
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto coarse_offset = coarse_cell.local_id * num_groups_;
    for (unsigned int g = 0; g < num_groups_; ++g)
      OpenSnPETScCall(VecSetValue(
        x, MapDOF(coarse_cell.global_id, g), coarse_phi_new_[coarse_offset + g], INSERT_VALUES));
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

  std::vector<int64_t> indices;
  indices.reserve(coarse_phi_new_.size());
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
    for (unsigned int g = 0; g < num_groups_; ++g)
      indices.push_back(MapDOF(coarse_cell.global_id, g));

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

  OpenSnInvalidArgumentIf(coarse_source_phi.size() != coarse_mesh_.NumLocalCells() * num_groups_,
                          "Coarse source scalar flux vector size mismatch.");

  OpenSnPETScCall(VecSet(rhs, 0.0));
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = *xs_map.at(coarse_cell.block_id);
    if (not xs.IsFissionable())
      continue;

    const auto& F = xs.GetProductionMatrix();
    const auto coarse_offset = coarse_cell.local_id * num_groups_;
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      double rhs_value = 0.0;
      for (unsigned int gp = 0; gp < num_groups_; ++gp)
      {
        rhs_value += F[g][gp] * coarse_source_phi[coarse_offset + gp];
      }
      rhs_value *= coarse_cell.volume / k_eff;
      const auto row = MapDOF(coarse_cell.global_id, g);
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
  OpenSnInvalidArgumentIf(coarse_phi.size() != coarse_mesh_.NumLocalCells() * num_groups_,
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
    const auto coarse_offset = coarse_cell.local_id * num_groups_;
    for (unsigned int g = 0; g < num_groups_; ++g)
      OpenSnPETScCall(VecSetValue(
        x, MapDOF(coarse_cell.global_id, g), coarse_phi[coarse_offset + g], INSERT_VALUES));
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

double
CMFDAcceleration::ComputeCoarseFissionProduction(const std::vector<double>& coarse_phi) const
{
  const auto& xs_map = do_problem_.GetBlockID2XSMap();

  OpenSnInvalidArgumentIf(coarse_phi.size() != coarse_mesh_.NumLocalCells() * num_groups_,
                          "Coarse scalar flux vector size mismatch.");

  double local_production = 0.0;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = *xs_map.at(coarse_cell.block_id);
    if (not xs.IsFissionable())
      continue;

    const auto& F = xs.GetProductionMatrix();
    const auto coarse_offset = coarse_cell.local_id * num_groups_;
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      for (unsigned int gp = 0; gp < num_groups_; ++gp)
      {
        local_production += F[g][gp] * coarse_phi[coarse_offset + gp] * coarse_cell.volume;
      }
    }
  }

  double global_production = 0.0;
  mpi_comm.all_reduce(local_production, global_production, mpi::op::sum<double>());
  return global_production;
}

void
CMFDAcceleration::LogTimingSummary(const double transport_time,
                                   const double restrict_old_time,
                                   const double restrict_new_time,
                                   const double coarse_flux_cache_time,
                                   const double face_current_cache_time,
                                   const double assemble_time,
                                   const double residual_time,
                                   const double fission_time,
                                   const double coarse_pi_time,
                                   const double update_time,
                                   const double post_time,
                                   const int coarse_pi_iterations,
                                   const int coarse_ksp_iterations,
                                   const BalanceResidual& residual) const
{
  const auto log_double_metric =
    [this](const std::string& category, const std::string& metric, const double value)
  {
    std::ostringstream out;
    out << std::scientific << std::setprecision(6);
    out << "CMFD_METRIC c=" << category << " i=" << outer_iteration_ << " m=" << metric
        << " v=" << value;
    log.Log() << program_timer.GetTimeString() << " " << out.str();
  };
  const auto log_int_metric =
    [this](const std::string& category, const std::string& metric, const int value)
  {
    log.Log() << program_timer.GetTimeString() << " CMFD_METRIC c=" << category
              << " i=" << outer_iteration_ << " m=" << metric << " v=" << value;
  };

  log_double_metric("timing", "transport_time", MaxAcrossRanks(transport_time));
  log_double_metric("timing", "restrict_old_time", MaxAcrossRanks(restrict_old_time));
  log_double_metric("timing", "restrict_new_time", MaxAcrossRanks(restrict_new_time));
  log_double_metric("timing", "coarse_flux_cache_time", MaxAcrossRanks(coarse_flux_cache_time));
  log_double_metric("timing", "face_current_cache_time", MaxAcrossRanks(face_current_cache_time));
  log_double_metric("timing", "assemble_time", MaxAcrossRanks(assemble_time));
  log_double_metric("timing", "residual_time", MaxAcrossRanks(residual_time));
  log_double_metric("timing", "fission_time", MaxAcrossRanks(fission_time));
  log_double_metric("timing", "coarse_pi_time", MaxAcrossRanks(coarse_pi_time));
  log_double_metric("timing", "update_time", MaxAcrossRanks(update_time));
  log_double_metric("timing", "post_time", MaxAcrossRanks(post_time));
  log_int_metric("timing", "coarse_pi_iterations", coarse_pi_iterations);
  log_int_metric("timing", "coarse_ksp_iterations", coarse_ksp_iterations);
  log_double_metric("timing", "balance_residual", residual.relative_l2);
}

CMFDAcceleration::CoarseMeshDiagnostics
CMFDAcceleration::ComputeCoarseMeshDiagnostics() const
{
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
  return diagnostics;
}

void
CMFDAcceleration::LogCoarseMeshDiagnostics(const CoarseMeshDiagnostics& diagnostics) const
{
  const double aggregation_ratio = diagnostics.global_coarse_cells > 0
                                     ? static_cast<double>(diagnostics.global_fine_cells) /
                                         static_cast<double>(diagnostics.global_coarse_cells)
                                     : 0.0;

  log.Log() << "CMFD_METRIC c=coarse_mesh m=global_fine_cells v=" << diagnostics.global_fine_cells;
  log.Log() << "CMFD_METRIC c=coarse_mesh m=global_coarse_cells v="
            << diagnostics.global_coarse_cells;
  log.Log() << "CMFD_METRIC c=coarse_mesh m=aggregation_ratio v=" << aggregation_ratio;
  log.Log() << "CMFD_METRIC c=coarse_mesh m=global_unknowns v=" << GlobalUnknownCount();
  log.Log() << "CMFD_METRIC c=coarse_mesh m=max_fine_cells_per_coarse_cell v="
            << diagnostics.global_max_fine_cells_per_coarse_cell;
  log.Log() << "CMFD_METRIC c=coarse_mesh m=max_faces_per_coarse_cell v="
            << diagnostics.global_max_faces_per_coarse_cell;
  log.Log() << "CMFD_METRIC c=coarse_mesh m=local_coarse_cells v="
            << diagnostics.local_coarse_cells;
  log.Log() << "CMFD_METRIC c=coarse_mesh m=undersized_coarse_cells v="
            << diagnostics.global_undersized_coarse_cells;
  const double average_faces_per_coarse_cell =
    diagnostics.global_coarse_cells > 0 ? diagnostics.global_faces_per_coarse_cell_sum /
                                            static_cast<double>(diagnostics.global_coarse_cells)
                                        : 0.0;
  log.Log() << "CMFD_METRIC c=coarse_mesh m=average_faces_per_coarse_cell v="
            << average_faces_per_coarse_cell;
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

bool
CMFDAcceleration::IsAcceptableFluxUpdate(const FluxUpdateDiagnostics& diagnostics,
                                         const double min_allowed_scalar_flux) const
{
  return not diagnostics.has_invalid_k and not diagnostics.has_nonfinite_flux and
         diagnostics.min_scalar_flux >= min_allowed_scalar_flux;
}

void
CMFDAcceleration::UpdateAdaptiveRelaxation(const double starting_relaxation,
                                           const double applied_damping,
                                           const bool skipped_correction)
{
  if (not adaptive_relaxation_)
    return;

  if (skipped_correction or
      applied_damping < adaptive_relaxation_accept_fraction_ * starting_relaxation)
  {
    accepted_strong_corrections_ = 0;
    const double reduced_relaxation = starting_relaxation * adaptive_relaxation_reduction_;
    current_relaxation_ =
      std::clamp(reduced_relaxation, adaptive_relaxation_min_, adaptive_relaxation_max_);
    return;
  }

  ++accepted_strong_corrections_;
  if (accepted_strong_corrections_ >= adaptive_relaxation_successes_to_grow_)
  {
    current_relaxation_ = std::clamp(starting_relaxation * adaptive_relaxation_growth_,
                                     adaptive_relaxation_min_,
                                     adaptive_relaxation_max_);
    accepted_strong_corrections_ = 0;
  }
}

double
CMFDAcceleration::ApplyFluxCorrectionWithDamping(const std::vector<double>& transport_phi_new,
                                                 const std::vector<double>& coarse_phi,
                                                 const double transport_k_eff,
                                                 const double cmfd_k_eff,
                                                 double& applied_damping,
                                                 unsigned int& attempts,
                                                 FluxUpdateDiagnostics& diagnostics,
                                                 bool& skipped_correction)
{
  std::vector<double> unrelaxed_coarse_delta_phi(coarse_phi.size(), 0.0);
  for (size_t i = 0; i < coarse_phi.size(); ++i)
    unrelaxed_coarse_delta_phi[i] = coarse_phi[i] - coarse_phi_new_[i];

  double damping = adaptive_relaxation_ ? current_relaxation_ : relaxation_;
  attempts = 0;
  skipped_correction = false;
  applied_damping = 0.0;
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

    auto candidate_phi_new = transport_phi_new;
    CMFDProlongateScalarFluxCorrection(
      do_problem_, first_group_, num_groups_, coarse_mesh_, coarse_delta_phi, candidate_phi_new);
    const double candidate_k_eff = transport_k_eff + damping * (cmfd_k_eff - transport_k_eff);
    diagnostics = AnalyzeFluxUpdate(candidate_phi_new, candidate_k_eff);
    if (IsAcceptableFluxUpdate(diagnostics, min_allowed_scalar_flux))
    {
      phi_new_local_ = std::move(candidate_phi_new);
      applied_damping = damping;
      ++attempts;
      return candidate_k_eff;
    }

    damping *= 0.5;
  }

  diagnostics = transport_diagnostics;
  phi_new_local_ = transport_phi_new;
  skipped_correction = true;
  if (verbose_ or outer_iteration_ < 5 or outer_iteration_ % 25 == 0)
    log.Log0Warning() << "CMFD correction rejected after " << attempts
                      << " damping attempts; using the unaccelerated transport update. "
                      << "Transport minimum scalar flux = " << diagnostics.min_scalar_flux << ".";
  return transport_k_eff;
}

void
CMFDAcceleration::LogCorrectionDiagnostics(const double starting_relaxation,
                                           const double applied_damping,
                                           const unsigned int attempts,
                                           const FluxUpdateDiagnostics& diagnostics,
                                           const bool skipped_correction) const
{
  log.Log() << "CMFD_METRIC c=correction i=" << outer_iteration_
            << " m=starting_relaxation v=" << starting_relaxation;
  log.Log() << "CMFD_METRIC c=correction i=" << outer_iteration_
            << " m=next_relaxation v=" << current_relaxation_;

  std::ostringstream out;
  out << std::scientific << std::setprecision(6);
  out << "CMFD_METRIC c=correction"
      << " i=" << outer_iteration_ << " m=damping v=" << applied_damping;
  log.Log() << out.str();

  log.Log() << "CMFD_METRIC c=correction i=" << outer_iteration_ << " m=attempts v=" << attempts;
  log.Log() << "CMFD_METRIC c=correction i=" << outer_iteration_
            << " m=skipped v=" << (skipped_correction ? 1 : 0);
  std::ostringstream min_flux;
  min_flux << std::scientific << std::setprecision(6);
  min_flux << "CMFD_METRIC c=correction"
           << " i=" << outer_iteration_ << " m=min_scalar_flux v=" << diagnostics.min_scalar_flux;
  log.Log() << min_flux.str();
  log.Log() << "CMFD_METRIC c=correction i=" << outer_iteration_
            << " m=nonfinite_flux v=" << (diagnostics.has_nonfinite_flux ? 1 : 0);
  log.Log() << "CMFD_METRIC c=correction i=" << outer_iteration_
            << " m=invalid_k v=" << (diagnostics.has_invalid_k ? 1 : 0);
}

} // namespace opensn
