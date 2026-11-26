// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_vector_tools.h"
#include "framework/parameters/input_parameters.h"
#include "framework/logging/log.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"
#include "framework/utils/timer.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/convergence.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace opensn
{
namespace
{

/// Consecutive skipped corrections after which CMFD warns that it is not accelerating.
constexpr unsigned int persistent_skip_warning_count = 5;
// A coarse-group flux sum within this fraction of the sum of its magnitudes is indistinguishable
// from roundoff, so homogenization falls back to volume weights for that coarse group.
constexpr double homogenization_cancellation_tolerance =
  64.0 * std::numeric_limits<double>::epsilon();
// A coarse linear solve is accepted if its true residual is at most this factor times the
// requested tolerance, max(coarse_linear_rel_tol * |b|, coarse_linear_abs_tol).
constexpr double true_residual_factor = 100.0;
// Tolerances and iteration limit of the coarse linear solves.
constexpr double coarse_linear_rel_tol = 1.0e-10;
constexpr double coarse_linear_abs_tol = 1.0e-200;
constexpr int coarse_linear_max_its = 1000;
// The coarse eigensolve: Wielandt-shifted power iteration whose flux and k tolerances are this
// fraction of the change the transport sweep made to the coarse flux, capped at the fixed
// tolerances below, with at most coarse_pi_max_its iterations.
constexpr double wielandt_shift = 0.05;
constexpr double coarse_pi_relative_tol = 1.0e-2;
constexpr double coarse_pi_phi_tol = 1.0e-6;
constexpr double coarse_pi_k_tol = 1.0e-8;
constexpr int coarse_pi_max_its = 1000;
// A coarse flux below this fraction of the largest flux of its coarse group counts as zero.
constexpr double zero_flux_fraction = 1.0e-14;
// Within-group tolerance of a transport update. Small enough that each update runs its full
// number of sweeps.
constexpr double transport_update_abs_tol = 1.0e-12;
// Number of coarse groups the automatic coarse_groups choice targets.
constexpr unsigned int auto_coarse_groups = 16;
// An LU factorization of the coarse operator is used while its estimated memory is at most this
// many bytes per MPI rank.
constexpr double coarse_direct_bytes_per_rank = 1.0 * 1024 * 1024 * 1024;
// Coarse cells thicker than this many transport mean free paths in some direction make CMFD with
// one sweep per update converge slowly.
constexpr double optically_thick_mfp = 3.0;

// Collective check: throws on every rank if any rank reports a positive count for a failure.
void
ThrowIfAnyRankCounted(const std::vector<std::pair<std::size_t, std::string>>& local_failures)
{
  std::vector<std::size_t> local_counts;
  local_counts.reserve(local_failures.size());
  for (const auto& failure : local_failures)
    local_counts.push_back(failure.first);
  std::vector<std::size_t> counts(local_counts.size(), 0);
  opensn::mpi_comm.all_reduce(local_counts, counts, mpi::op::sum<std::size_t>());
  for (std::size_t i = 0; i < counts.size(); ++i)
    if (counts[i] > 0)
      throw std::logic_error("CMFD homogenization: " + std::to_string(counts[i]) + " " +
                             local_failures[i].second + ".");
}

// Whether an MPI all-to-all buffer of this many values fits the int counts MPI uses.
bool
FitsMPICount(const std::size_t count)
{
  return count <= static_cast<std::size_t>(std::numeric_limits<int>::max());
}

std::vector<double>
BuildUniformEdges(const double origin, const double pitch, const int count)
{
  std::vector<double> edges(static_cast<std::size_t>(count) + 1);
  for (int i = 0; i <= count; ++i)
    edges[static_cast<std::size_t>(i)] = origin + static_cast<double>(i) * pitch;
  return edges;
}

/**
 * Half transmissibility of a coarse cell to one of its faces, D A / |d|, with d the vector from the
 * cell centroid to the face centroid. It is the usual D A / (d.n) for a face the centroid projects
 * onto orthogonally. On irregular aggregates d can be nearly parallel to the face, where d.n would
 * make the coupling grow without bound; |d| keeps it bounded without weakening it on skewed faces.
 *
 * This is the orthogonal-correction decomposition of the face-normal gradient used by finite-volume
 * methods on non-orthogonal meshes (Jasak, PhD thesis, Imperial College, 1996; Moukalled, Mangani
 * and Darwish, The Finite Volume Method in Computational Fluid Dynamics, Springer, 2016, ch. 8).
 * The non-orthogonal remainder, which those methods add as an explicit correction, needs no
 * separate term here: the CMFD current correction (D-hat) makes the coarse currents reproduce the
 * transport currents for any coupling, so the coupling affects only the stability of the iteration.
 */
double
HalfTransmissibility(const double diffusion_coefficient,
                     const Vector3& cell_centroid,
                     const Vector3& face_centroid,
                     const double face_area)
{
  const double distance = (face_centroid - cell_centroid).Norm();
  OpenSnLogicalErrorIf(distance <= 0.0,
                       "CMFD coarse face centroid coincides with its cell centroid.");
  return diffusion_coefficient * face_area / distance;
}

/// Diffusion coupling of two coarse cells across a face: the harmonic combination of the two
/// half transmissibilities.
double
InterfaceDiffusionCoefficient(const double owner_half_transmissibility,
                              const double neighbor_half_transmissibility)
{
  const double sum = owner_half_transmissibility + neighbor_half_transmissibility;
  return sum > 0.0 ? owner_half_transmissibility * neighbor_half_transmissibility / sum : 0.0;
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

/**
 * Relative L2 change between two coarse scalar-flux iterates, each normalized by its fission
 * production so that the arbitrary eigenvector scale of power iteration does not enter.
 */
double
NormalizedCoarseFluxChange(const std::vector<double>& phi_old,
                           const double production_old,
                           const std::vector<double>& phi_new,
                           const double production_new)
{
  double local_sums[2] = {0.0, 0.0};
  for (std::size_t i = 0; i < phi_new.size(); ++i)
  {
    const double shape_new = phi_new[i] / production_new;
    const double difference = shape_new - phi_old[i] / production_old;
    local_sums[0] += difference * difference;
    local_sums[1] += shape_new * shape_new;
  }
  double global_sums[2] = {0.0, 0.0};
  opensn::mpi_comm.all_reduce(local_sums, 2, global_sums, mpi::op::sum<double>());
  return std::sqrt(global_sums[0] / std::max(global_sums[1], 1.0e-300));
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

InputParameters
CMFDAcceleration::GetInputParameters()
{
  auto params = DiscreteOrdinatesKEigenAcceleration::GetInputParameters();

  params.AddOptionalParameter(
    "coarse_mesh",
    std::string("local_aggregation"),
    "CMFD coarse mesh. A string selects a type with its defaults, and a dict selects a type with "
    "its parameters. The type \"local_aggregation\" (dict key \"size\", the target number of "
    "fine cells per coarse cell, default 32) builds aggregates within each MPI rank. The type "
    "\"global_aggregation\" (key \"size\") builds connected aggregates that may span ranks. The "
    "type \"reactor_lattice\" (keys \"origin\" [x, y], \"pitch\" [px, py], \"shape\" [nx, ny] "
    "and optional \"axial_edges\") gives one coarse cell per lattice cell and axial layer. The "
    "type \"cartesian\" (keys \"x_edges\", \"y_edges\", \"z_edges\", each optional but not all "
    "empty) gives one coarse cell per Cartesian box. The type \"logical_volumes\" (key "
    "\"volumes\") gives one coarse cell per connected part of each volume. The type "
    "\"identity\" gives one coarse cell per transport cell and is meant for testing. Coarse "
    "cells may contain several materials. The coarse-cell coefficients are homogenized from "
    "the transport flux.");
  params.SetParameterTypeMismatchAllowed("coarse_mesh");
  params.AddOptionalParameter(
    "coarse_groups",
    std::string("auto"),
    "Number of CMFD energy groups, or \"auto\". The transport groups are split into consecutive "
    "coarse groups of ceil(G / coarse_groups) fine groups each, which may give fewer coarse "
    "groups than requested. The choice \"auto\" keeps every group when G <= 16 and targets 16 "
    "coarse groups otherwise.");
  params.SetParameterTypeMismatchAllowed("coarse_groups");

  // Names released before the CMFD interface was simplified. Each is accepted for one release
  // with a warning that names the replacement. Values that contradict the fixed method are errors
  // in ParseCoarseMeshSpec or CheckRequirements.
  params.AddOptionalParameter("aggregation_size", 32, "Deprecated. Use coarse_mesh.");
  params.MarkParameterDeprecatedWarning(
    "aggregation_size", R"(Use coarse_mesh={"type": "local_aggregation", "size": n}.)");
  params.AddOptionalParameter("group_aggregation_size", 1, "Deprecated. Use coarse_groups.");
  params.MarkParameterDeprecatedWarning("group_aggregation_size",
                                        "Use coarse_groups = ceil(G / group_aggregation_size).");
  params.AddOptionalParameter(
    "current_closure", std::string("partial"), "Deprecated. CMFD always uses the partial closure.");
  params.MarkParameterDeprecatedWarning(
    "current_closure", "CMFD always uses the partial-current closure, and \"net\" is an error.");
  params.AddOptionalParameter("relaxation", 1.0, "Deprecated. CMFD corrections are unrelaxed.");
  params.MarkParameterDeprecatedWarning("relaxation",
                                        "CMFD corrections are unrelaxed, and a value other "
                                        "than 1 is an error.");
  params.AddOptionalParameter("correction_max_attempts", 10, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning(
    "correction_max_attempts",
    "CMFD no longer damps corrections. A correction that is not admissible is skipped.");
  params.AddOptionalParameter("correction_min_damping", 1.0e-4, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning(
    "correction_min_damping",
    "CMFD no longer damps corrections. A correction that is not admissible is skipped.");
  params.AddOptionalParameter("inactive_iterations", 0, "Deprecated. CMFD starts at once.");
  params.MarkParameterDeprecatedWarning(
    "inactive_iterations", "CMFD corrects every update, and a value other than 0 is an error.");
  params.AddOptionalParameter(
    "update_wgs_max_its", 1, "Deprecated. Use sweeps_per_update of the power iteration solver.");
  params.MarkParameterDeprecatedWarning("update_wgs_max_its",
                                        "Use sweeps_per_update of PowerIterationKEigenSolver.");
  params.AddOptionalParameter("update_wgs_abs_tol", 1.0e-12, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning("update_wgs_abs_tol",
                                        "Use sweeps_per_update of PowerIterationKEigenSolver.");
  params.AddOptionalParameter("balance_residual_tolerance", 1.0e-6, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning(
    "balance_residual_tolerance",
    "Use fission_source_tol of PowerIterationKEigenSolver, which tests the converged shape.");
  params.AddOptionalParameter("coarse_solver_policy", std::string("auto"), "Deprecated.");
  params.MarkParameterDeprecatedWarning("coarse_solver_policy",
                                        "CMFD chooses the coarse linear solver. Use "
                                        "petsc_options to configure it instead.");
  params.AddOptionalParameter("coarse_direct_solve_threshold", 150000, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning(
    "coarse_direct_solve_threshold",
    "CMFD chooses the coarse linear solver from its estimated memory.");
  params.AddOptionalParameter("l_rel_tol", 1.0e-8, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning("l_rel_tol",
                                        "CMFD sets the coarse linear solver tolerances.");
  params.ChangeExistingParamToOptional("l_abs_tol", 1.0e-7, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning("l_abs_tol",
                                        "CMFD sets the coarse linear solver tolerances.");
  params.ChangeExistingParamToOptional("max_iters", 100, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning("max_iters",
                                        "CMFD sets the coarse linear solver iteration limit.");
  params.ChangeExistingParamToOptional("pi_max_its", 1000, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning("pi_max_its",
                                        "CMFD converges the coarse eigenproblem itself.");
  params.ChangeExistingParamToOptional("pi_k_tol", 1.0e-8, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning("pi_k_tol",
                                        "CMFD converges the coarse eigenproblem itself.");
  params.AddOptionalParameter("pi_phi_tol", 1.0e-6, "Deprecated and ignored.");
  params.MarkParameterDeprecatedWarning("pi_phi_tol",
                                        "CMFD converges the coarse eigenproblem itself.");

  params.ConstrainParameterRange("aggregation_size", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("group_aggregation_size", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("current_closure", AllowableRangeList::New({"net", "partial"}));
  params.ConstrainParameterRange("update_wgs_max_its", AllowableRangeLowLimit::New(1));
  params.ChangeExistingParamToOptional("name", "CMFDAcceleration");
  params.ChangeExistingParamToOptional(
    "verbose", false, "If true, prints CMFD coarse-mesh, correction and timing diagnostics.");
  params.ChangeExistingParamToOptional(
    "petsc_options",
    std::string(""),
    "PETSc options for the coarse linear solver, prefixed by the accelerator name. Options "
    "that set the KSP or PC type replace the automatic choice of the coarse linear solver.");

  return params;
}

std::shared_ptr<CMFDAcceleration>
CMFDAcceleration::Create(const ParameterBlock& params)
{
  return CreateObject<CMFDAcceleration>("lbs::CMFDAcceleration", params);
}

namespace
{

bool
WasGiven(const InputParameters& params, const std::string& name)
{
  return params.GetParametersAtAssignment().Has(name);
}

std::vector<double>
OptionalDoubles(const ParameterBlock& block, const std::string& key)
{
  return block.Has(key) ? block.GetParamVectorValue<double>(key) : std::vector<double>{};
}

} // namespace

CMFDAcceleration::CoarseMeshSpec
CMFDAcceleration::ParseCoarseMeshSpec(const InputParameters& params)
{
  CoarseMeshSpec spec;
  const auto& coarse_mesh = params.GetParam("coarse_mesh");
  const auto Fail = [](const std::string& message)
  { throw std::invalid_argument("CMFDAcceleration coarse_mesh: " + message); };

  if (coarse_mesh.GetType() == ParameterBlockType::STRING)
    spec.type = coarse_mesh.GetValue().GetStringValue();
  else if (coarse_mesh.GetType() == ParameterBlockType::BLOCK)
  {
    if (not coarse_mesh.Has("type"))
      Fail("a dict needs a \"type\" key.");
    spec.type = coarse_mesh.GetParamValue<std::string>("type");
  }
  else
    Fail("expected a string or a dict.");

  const std::map<std::string, std::set<std::string>> allowed_keys = {
    {"identity", {"type"}},
    {"local_aggregation", {"type", "size"}},
    {"global_aggregation", {"type", "size"}},
    {"reactor_lattice", {"type", "origin", "pitch", "shape", "axial_edges"}},
    {"cartesian", {"type", "x_edges", "y_edges", "z_edges"}},
    {"logical_volumes", {"type", "volumes"}}};
  const auto allowed = allowed_keys.find(spec.type);
  if (allowed == allowed_keys.end())
    Fail("unknown type \"" + spec.type +
         "\". Valid types are identity, local_aggregation, global_aggregation, "
         "reactor_lattice, cartesian, and logical_volumes.");

  if (coarse_mesh.GetType() == ParameterBlockType::BLOCK)
  {
    for (const auto& entry : coarse_mesh.GetParameters())
      if (allowed->second.count(entry.GetName()) == 0)
        Fail("type \"" + spec.type + "\" does not take the key \"" + entry.GetName() + "\".");
    if (coarse_mesh.Has("size"))
      spec.size = coarse_mesh.GetParamValue<int>("size");
    spec.x_edges = OptionalDoubles(coarse_mesh, "x_edges");
    spec.y_edges = OptionalDoubles(coarse_mesh, "y_edges");
    spec.z_edges = OptionalDoubles(coarse_mesh, "z_edges");
    spec.lattice_origin = OptionalDoubles(coarse_mesh, "origin");
    spec.lattice_pitch = OptionalDoubles(coarse_mesh, "pitch");
    spec.lattice_axial_edges = OptionalDoubles(coarse_mesh, "axial_edges");
    if (coarse_mesh.Has("shape"))
      spec.lattice_shape = coarse_mesh.GetParamVectorValue<int>("shape");
    if (coarse_mesh.Has("volumes"))
      spec.logical_volumes =
        coarse_mesh.GetParamVectorValue<std::shared_ptr<LogicalVolume>>("volumes");
    if (spec.type == "reactor_lattice" and
        (spec.lattice_origin.empty() or spec.lattice_pitch.empty() or spec.lattice_shape.empty()))
      Fail(R"(type "reactor_lattice" needs "origin", "pitch", and "shape".)");
    if (spec.type == "logical_volumes" and spec.logical_volumes.empty())
      Fail(R"(type "logical_volumes" needs a nonempty "volumes" list.)");
    if (spec.type == "cartesian" and spec.x_edges.empty() and spec.y_edges.empty() and
        spec.z_edges.empty())
      Fail(R"(type "cartesian" needs at least one of "x_edges", "y_edges", "z_edges".)");
  }
  else if (spec.type == "reactor_lattice" or spec.type == "cartesian" or
           spec.type == "logical_volumes")
    Fail("type \"" + spec.type + "\" needs a dict with its geometry.");

  if (WasGiven(params, "aggregation_size"))
  {
    if (spec.type != "local_aggregation" and spec.type != "global_aggregation")
      Fail("the deprecated aggregation_size applies only to the aggregation types.");
    if (coarse_mesh.GetType() == ParameterBlockType::BLOCK and coarse_mesh.Has("size"))
      Fail("give the aggregation size once, as the \"size\" key.");
    spec.size = params.GetParamValue<int>("aggregation_size");
  }
  if (spec.size < 1)
    Fail("the aggregation \"size\" must be at least 1.");
  return spec;
}

unsigned int
CMFDAcceleration::ParseGroupAggregationSize(const InputParameters& params,
                                            const unsigned int num_groups)
{
  const auto& coarse_groups = params.GetParam("coarse_groups");
  const bool given = WasGiven(params, "coarse_groups");
  if (WasGiven(params, "group_aggregation_size"))
  {
    if (given)
      throw std::invalid_argument("CMFDAcceleration: give coarse_groups or the deprecated "
                                  "group_aggregation_size, not both.");
    return static_cast<unsigned int>(params.GetParamValue<int>("group_aggregation_size"));
  }
  unsigned int target = 0;
  if (coarse_groups.GetType() == ParameterBlockType::STRING)
  {
    if (coarse_groups.GetValue().GetStringValue() != "auto")
      throw std::invalid_argument(
        "CMFDAcceleration coarse_groups: expected a positive integer or \"auto\".");
    target = num_groups <= auto_coarse_groups ? num_groups : auto_coarse_groups;
  }
  else if (coarse_groups.GetType() == ParameterBlockType::INTEGER)
  {
    const auto value = coarse_groups.GetValue().GetIntegerValue();
    if (value < 1)
      throw std::invalid_argument("CMFDAcceleration coarse_groups: must be at least 1.");
    target = static_cast<unsigned int>(std::min<int64_t>(value, num_groups));
  }
  else
    throw std::invalid_argument(
      "CMFDAcceleration coarse_groups: expected a positive integer or \"auto\".");
  return (num_groups + target - 1) / target;
}

CMFDAcceleration::CMFDAcceleration(const InputParameters& params)
  : DiscreteOrdinatesKEigenAcceleration(params),
    coarse_mesh_spec_(ParseCoarseMeshSpec(params)),
    legacy_update_wgs_max_its_(
      WasGiven(params, "update_wgs_max_its")
        ? static_cast<unsigned int>(params.GetParamValue<int>("update_wgs_max_its"))
        : 0),
    num_groups_(do_problem_.GetNumGroups()),
    group_aggregation_size_(ParseGroupAggregationSize(params, num_groups_)),
    num_coarse_groups_((num_groups_ + group_aggregation_size_ - 1) / group_aggregation_size_)
{
  const auto Reject = [](const std::string& message)
  { throw std::invalid_argument("CMFDAcceleration: " + message); };
  if (params.GetParamValue<std::string>("current_closure") != "partial")
    Reject("current_closure=\"net\" is no longer supported. CMFD uses the partial-current "
           "closure.");
  if (params.GetParamValue<double>("relaxation") != 1.0)
    Reject("relaxation is fixed at 1. Remove the parameter.");
  if (params.GetParamValue<int>("inactive_iterations") != 0)
    Reject("inactive_iterations is no longer supported. Remove the parameter.");
}

CMFDAcceleration::~CMFDAcceleration()
{
  if (ksp_)
    OpenSnPETScCall(KSPDestroy(&ksp_));
  if (rhs_)
    OpenSnPETScCall(VecDestroy(&rhs_));
  if (A_)
    OpenSnPETScCall(MatDestroy(&A_));
  if (A_shift_)
    OpenSnPETScCall(MatDestroy(&A_shift_));
}

const CMFDAcceleration::CoarseCellXS&
CMFDAcceleration::LocalHomogenizedXS(const uint64_t coarse_cell_global_id) const
{
  OpenSnLogicalErrorIf(coarse_mesh_.NumGlobalCells() == 0 or
                         homogenized_xs_.size() != coarse_mesh_.NumLocalCells(),
                       "CMFD homogenized coefficients are not available before the owning "
                       "solver's first Execute().");
  OpenSnInvalidArgumentIf(not coarse_mesh_.HasLocalCoarseCell(coarse_cell_global_id),
                          "CMFD coarse cell " + std::to_string(coarse_cell_global_id) +
                            " is not owned by rank " + std::to_string(opensn::mpi_comm.rank()) +
                            ".");
  return homogenized_xs_[coarse_mesh_.LocalCellFromGlobalID(coarse_cell_global_id).local_id];
}

std::vector<double>
CMFDAcceleration::GetHomogenizedRemovalXS(const uint64_t coarse_cell_global_id) const
{
  return LocalHomogenizedXS(coarse_cell_global_id).removal;
}

std::vector<double>
CMFDAcceleration::GetHomogenizedDiffusionCoefficient(const uint64_t coarse_cell_global_id) const
{
  return LocalHomogenizedXS(coarse_cell_global_id).diffusion;
}

std::vector<double>
CMFDAcceleration::GetHomogenizedScatteringXS(const uint64_t coarse_cell_global_id) const
{
  std::vector<double> matrix;
  ComputeCoarseScattering(LocalHomogenizedXS(coarse_cell_global_id), matrix);
  return matrix;
}

std::vector<double>
CMFDAcceleration::GetHomogenizedProductionXS(const uint64_t coarse_cell_global_id) const
{
  std::vector<double> matrix;
  ComputeCoarseProduction(LocalHomogenizedXS(coarse_cell_global_id), matrix);
  return matrix;
}

void
CMFDAcceleration::CheckRequirements(std::vector<std::string>& errors) const
{
  // The coarse balance uses mesh cell volumes and face areas. In 1D spherical geometry the face
  // areas are spherical (zero at r = 0) but the volumes are slab lengths, so the coarse system is
  // inconsistent: aggregated coarse meshes have a degenerate face and identity meshes do not
  // accelerate. Homogenized coefficients are ratios of integrals in one measure, so they do not
  // add a geometry restriction.
  if (do_problem_.GetGeometryType() == GeometryType::ONED_SPHERICAL)
    errors.emplace_back(GetName() + " is not supported in 1D spherical geometry.");

  const auto& spec = coarse_mesh_spec_;
  const auto CheckEdges = [&](const std::vector<double>& edges, const std::string& name)
  {
    if (edges.empty())
      return;
    if (edges.size() < 2)
    {
      errors.emplace_back(GetName() + " requires at least two " + name + " values.");
      return;
    }
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
      if (not std::isfinite(edges[i]))
      {
        errors.emplace_back(GetName() + " requires finite " + name + " values.");
        return;
      }
      if (i > 0 and edges[i] <= edges[i - 1])
      {
        errors.emplace_back(GetName() + " requires strictly increasing " + name + ".");
        return;
      }
    }
  };
  if (spec.type == "cartesian")
  {
    CheckEdges(spec.x_edges, "x_edges");
    CheckEdges(spec.y_edges, "y_edges");
    CheckEdges(spec.z_edges, "z_edges");
  }
  else if (spec.type == "reactor_lattice")
  {
    if (spec.lattice_origin.size() != 2)
      errors.emplace_back(GetName() + " requires two reactor_lattice origin values (x and y).");
    if (spec.lattice_pitch.size() != 2 or
        std::any_of(spec.lattice_pitch.begin(),
                    spec.lattice_pitch.end(),
                    [](const double value) { return not std::isfinite(value) or value <= 0.0; }))
      errors.emplace_back(GetName() + " requires two finite positive reactor_lattice pitches.");
    if (spec.lattice_shape.size() != 2 or std::any_of(spec.lattice_shape.begin(),
                                                      spec.lattice_shape.end(),
                                                      [](const int value) { return value <= 0; }))
      errors.emplace_back(GetName() + " requires two positive reactor_lattice shape values.");
    if (std::any_of(spec.lattice_origin.begin(),
                    spec.lattice_origin.end(),
                    [](const double value) { return not std::isfinite(value); }))
      errors.emplace_back(GetName() + " requires finite reactor_lattice origin values.");
    if (spec.lattice_origin.size() == 2 and spec.lattice_pitch.size() == 2 and
        spec.lattice_shape.size() == 2)
      for (std::size_t d = 0; d < 2; ++d)
        if (not std::isfinite(spec.lattice_origin[d] +
                              spec.lattice_pitch[d] * static_cast<double>(spec.lattice_shape[d])))
          errors.emplace_back(GetName() + " reactor_lattice extent must be finite.");
    CheckEdges(spec.lattice_axial_edges, "reactor_lattice axial_edges");
  }

  const bool adjoint = do_problem_.GetOptions().adjoint;
  for (const auto& [block_id, xs] : do_problem_.GetBlockID2XSMap())
  {
    if (not xs->DiffusionInitialized() or xs->GetSigmaTransport(adjoint).size() < num_groups_)
    {
      errors.emplace_back(GetName() + " requires diffusion data for material block " +
                          std::to_string(block_id) + ".");
      continue;
    }
    for (unsigned int g = 0; g < num_groups_; ++g)
      if (const double sigma_tr = xs->GetSigmaTransport(adjoint)[g];
          not std::isfinite(sigma_tr) or sigma_tr < 0.0)
      {
        errors.emplace_back(GetName() +
                            " requires finite non-negative transport cross sections; "
                            "material block " +
                            std::to_string(block_id) + ", group " + std::to_string(g) +
                            " is invalid.");
        break;
      }
  }
}

void
CMFDAcceleration::Initialize()
{
  do_problem_.ConfigureOutflowStorage(true);
  // Coarse vectors from an earlier Execute() are indexed by the previous coarse mesh.
  coarse_phi_old_.clear();
  coarse_phi_new_.clear();

  // Every coarse-mesh type may combine materials in one coarse cell. Homogenization weights each
  // material by its transport flux, so that the coarse cell preserves the reaction rates.
  const auto& spec = coarse_mesh_spec_;
  const auto& grid = *do_problem_.GetGrid();
  if (spec.type == "identity")
    coarse_mesh_ = CMFDCoarseMesh::BuildIdentity(grid);
  else if (spec.type == "local_aggregation")
    coarse_mesh_ = CMFDCoarseMesh::BuildLocalAggregation(grid, spec.size, true);
  else if (spec.type == "global_aggregation")
    coarse_mesh_ = CMFDCoarseMesh::BuildGlobalAggregation(grid, spec.size, true);
  else if (spec.type == "logical_volumes")
    coarse_mesh_ = CMFDCoarseMesh::BuildLogicalVolumeOverlay(grid, spec.logical_volumes);
  else if (spec.type == "cartesian")
    coarse_mesh_ =
      CMFDCoarseMesh::BuildCartesianOverlay(grid, spec.x_edges, spec.y_edges, spec.z_edges);
  else if (spec.type == "reactor_lattice")
    coarse_mesh_ = CMFDCoarseMesh::BuildCartesianOverlay(
      grid,
      BuildUniformEdges(
        spec.lattice_origin.at(0), spec.lattice_pitch.at(0), spec.lattice_shape.at(0)),
      BuildUniformEdges(
        spec.lattice_origin.at(1), spec.lattice_pitch.at(1), spec.lattice_shape.at(1)),
      spec.lattice_axial_edges);
  else
    OpenSnLogicalError("Unknown CMFD coarse mesh type " + spec.type + ".");

  if (const auto& diagnostics = coarse_mesh_.OverlayDiagnostics())
  {
    const double area_fraction = diagnostics->interface_area > 0.0
                                   ? diagnostics->nonconforming_area / diagnostics->interface_area
                                   : 0.0;
    log.Log() << "CMFD overlay: " << coarse_mesh_.NumGlobalCells() << " coarse cells. "
              << diagnostics->num_nonconforming_faces << " of " << diagnostics->num_interface_faces
              << " coarse-interface fine faces do not lie on an overlay plane (area fraction "
              << area_fraction << ").";
    if (diagnostics->num_thin_intervals > 0)
      log.Log0Warning() << "CMFD overlay has " << diagnostics->num_thin_intervals
                        << " occupied intervals narrower than "
                        << CMFDCoarseMesh::overlay_sliver_fraction
                        << " of the median occupied width of their axis (smallest ratio "
                        << diagnostics->min_interval_width_ratio
                        << "). Thin coarse cells give poorly scaled coarse-operator rows.";
  }

  InitializeLinearSystem();
  InitializeHomogenization();
  BuildHomogenizedXS({}, false);
  BuildCoarseFluxCache();
  AssembleOperator(false);

  CoarseMeshDiagnostics diagnostics;
  diagnostics.local_coarse_cells = coarse_mesh_.NumLocalCells();

  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    diagnostics.local_fine_cells += coarse_cell.fine_cell_ids.size();
    diagnostics.local_faces_per_coarse_cell_sum += static_cast<double>(coarse_cell.faces.size());
    if (coarse_cell.fine_cell_ids.size() < static_cast<std::size_t>(spec.size))
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
CMFDAcceleration::PreExecute()
{
  coarse_phi_old_.clear();
  coarse_phi_new_.clear();
  coarse_phi_cache_.clear();
  face_current_cache_.clear();
  outer_iteration_ = 0;
  num_applied_corrections_ = 0;
  num_skipped_corrections_ = 0;
  consecutive_skipped_corrections_ = 0;
  warned_persistent_skips_ = false;
  warned_homogenization_fallback_ = false;
  warned_coarse_pi_limit_ = false;
  warned_kept_coarse_values_ = false;
  last_kept_coarse_values_ = 0;
  checked_optical_thickness_ = false;
  num_closure_fallbacks_ = 0;
  sweeps_per_update_ = legacy_update_wgs_max_its_ > 0 and solver_->IsSweepsPerUpdateAutomatic()
                         ? legacy_update_wgs_max_its_
                         : solver_->GetSweepsPerUpdate();
  last_restrict_old_time_ = 0.0;
  transport_start_time_ = 0.0;
  old_fission_production_ = 0.0;
}

void
CMFDAcceleration::PrePowerIteration()
{
  // The owning solver starts each Execute() with the caller's WGS settings and restores them
  // when Execute() returns. Every transport update between CMFD corrections is
  // sweeps_per_update_ within-group iterations.
  if (outer_iteration_ == 0)
  {
    ConfigureTransportSolve(sweeps_per_update_, transport_update_abs_tol);
    log.Log() << no_wrap << program_timer.GetTimeString() << " CMFD switching " << groupsets_.size()
              << " groupset WGS solver(s) to " << sweeps_per_update_
              << " iteration(s) per transport update.";
  }

  const double t0 = NowSeconds();
  // Fail fast if the pre-sweep flux is already degenerate, before paying for the sweep. This is
  // also the reference for an unaccelerated power update: WGS overwrites phi_old_local_.
  old_fission_production_ = ComputeFissionProduction(do_problem_, phi_old_local_);
  OpenSnLogicalErrorIf(old_fission_production_ == 0.0,
                       "CMFDAcceleration cannot update k with zero old production.");
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
  coarse_phi_new_ = CMFDRestrictScalarFlux(
    do_problem_, first_group_, num_groups_, group_aggregation_size_, coarse_mesh_, phi_new_local_);
  const double restrict_new_time = NowSeconds() - t0;

  t0 = NowSeconds();
  const double F_transport_new = ComputeFissionProduction(do_problem_, phi_new_local_);
  const double raw_transport_k_eff =
    F_transport_new / old_fission_production_ * solver_->GetEigenvalue();
  // Corrections use the post-WGS production for their normalization. Keep this separate from the
  // pre-sweep production used by a skipped correction.
  const double cmfd_reference_fission_production =
    ComputeFissionProduction(do_problem_, phi_old_local_);
  OpenSnLogicalErrorIf(cmfd_reference_fission_production == 0.0,
                       "CMFDAcceleration cannot update k with zero CMFD reference production.");
  const double transport_k_eff =
    F_transport_new / cmfd_reference_fission_production * solver_->GetEigenvalue();
  const double fission_time = NowSeconds() - t0;

  t0 = NowSeconds();
  BuildHomogenizedXS(phi_new_local_, true);
  const double homogenization_time = NowSeconds() - t0;

  t0 = NowSeconds();
  BuildCoarseFluxCache();
  const double coarse_flux_cache_time = NowSeconds() - t0;

  t0 = NowSeconds();
  BuildFaceCurrentCache();
  const double face_current_cache_time = NowSeconds() - t0;

  if (not checked_optical_thickness_)
  {
    CheckOpticalThickness();
    checked_optical_thickness_ = true;
  }

  t0 = NowSeconds();
  const auto sign_counts = AssembleOperator(true);
  const double assemble_time = NowSeconds() - t0;
  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
  {
    // Signed transport iterates can give the coarse operator entries of the wrong sign for an
    // M-matrix, which is what makes coarse solves unreliable. These counts help diagnose that.
    const std::array<std::size_t, 3> local_counts = {
      sign_counts.nonpositive_diagonals,
      sign_counts.positive_off_diagonals,
      static_cast<std::size_t>(std::count_if(
        coarse_phi_new_.begin(), coarse_phi_new_.end(), [](double phi) { return phi < 0.0; }))};
    std::array<std::size_t, 3> global_counts = {0, 0, 0};
    for (std::size_t i = 0; i < local_counts.size(); ++i)
      mpi_comm.all_reduce(local_counts[i], global_counts[i], mpi::op::sum<std::size_t>());
    log.Log() << no_wrap << program_timer.GetTimeString()
              << " CMFD operator signs: nonpositive_diagonals = " << global_counts[0]
              << ", positive_off_diagonals = " << global_counts[1]
              << ", negative_coarse_fluxes = " << global_counts[2]
              << ", rows = " << GlobalUnknownCount();
    // Negative restricted transport fluxes by coarse group, with the most negative value relative
    // to the group's largest flux.
    std::vector<double> local_group(2 * static_cast<std::size_t>(num_coarse_groups_), 0.0);
    for (std::size_t i = 0; i < coarse_phi_new_.size(); ++i)
    {
      const auto cg = i % num_coarse_groups_;
      if (coarse_phi_new_[i] < 0.0)
        local_group[cg] += 1.0;
      local_group[num_coarse_groups_ + cg] =
        std::max(local_group[num_coarse_groups_ + cg], std::fabs(coarse_phi_new_[i]));
    }
    std::vector<double> local_min(num_coarse_groups_, 0.0);
    for (std::size_t i = 0; i < coarse_phi_new_.size(); ++i)
      local_min[i % num_coarse_groups_] =
        std::min(local_min[i % num_coarse_groups_], coarse_phi_new_[i]);
    std::vector<double> global_min(num_coarse_groups_, 0.0);
    std::vector<double> counts_only(local_group.begin(), local_group.begin() + num_coarse_groups_);
    std::vector<double> max_only(local_group.begin() + num_coarse_groups_, local_group.end());
    std::vector<double> global_counts_g(num_coarse_groups_, 0.0);
    std::vector<double> global_max_g(num_coarse_groups_, 0.0);
    mpi_comm.all_reduce(counts_only.data(),
                        static_cast<int>(num_coarse_groups_),
                        global_counts_g.data(),
                        mpi::op::sum<double>());
    mpi_comm.all_reduce(max_only.data(),
                        static_cast<int>(num_coarse_groups_),
                        global_max_g.data(),
                        mpi::op::max<double>());
    mpi_comm.all_reduce(local_min.data(),
                        static_cast<int>(num_coarse_groups_),
                        global_min.data(),
                        mpi::op::min<double>());
    std::ostringstream by_group;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      by_group << " " << cg << ":" << static_cast<std::size_t>(global_counts_g[cg]) << "("
               << std::setprecision(2) << global_min[cg] / std::max(global_max_g[cg], 1.0e-300)
               << ")";
    log.Log() << no_wrap << program_timer.GetTimeString()
              << " CMFD negative coarse fluxes by coarse group (count and most negative relative "
                 "to the group maximum):"
              << by_group.str();
  }

  t0 = NowSeconds();
  BalanceResidual transport_balance_residual;
  BalanceResidual transport_current_balance_residual;
  if (verbose_)
  {
    transport_balance_residual =
      ComputeBalanceResidual(transport_k_eff, coarse_phi_new_, coarse_phi_old_);
    transport_current_balance_residual =
      ComputeTransportCurrentBalanceResidual(transport_k_eff, coarse_phi_new_, coarse_phi_old_);
    last_balance_residual_known_ = true;
    last_transport_current_balance_relative_l2_ = transport_current_balance_residual.relative_l2;
    if (outer_iteration_ < 5 or outer_iteration_ % 25 == 0)
    {
      log.Log() << no_wrap << program_timer.GetTimeString()
                << " CMFD operator balance residual: max_abs = "
                << transport_balance_residual.max_abs
                << ", relative_l2 = " << transport_balance_residual.relative_l2;
      log.Log() << no_wrap << program_timer.GetTimeString()
                << " CMFD transport-current balance residual: max_abs = "
                << transport_current_balance_residual.max_abs
                << ", relative_l2 = " << transport_current_balance_residual.relative_l2;
    }
  }
  const double residual_time = NowSeconds() - t0;

  auto coarse_phi_m = coarse_phi_new_;
  auto coarse_phi = coarse_phi_new_;
  double k_eff = transport_k_eff;
  double coarse_production_old = ComputeCoarseFissionProduction(coarse_phi_m);
  const double initial_coarse_production = coarse_production_old;
  // The coarse eigensolve breaks down when the coarse production or k of an iterate is not a
  // positive finite number. The correction is then skipped. The production is a global sum, so
  // every rank takes the same branch.
  bool coarse_eigensolve_failed =
    not std::isfinite(coarse_production_old) or coarse_production_old <= 0.0;
  // Tolerances of this coarse eigensolve: a fraction of the change the transport sweep made to the
  // production-normalized coarse flux, capped at fixed values. A fixed tolerance alone would add
  // an error of its size to every correction, so the outer iteration would stall once its own
  // error reached that level.
  const double transport_coarse_change =
    NormalizedCoarseFluxChange(coarse_phi_old_,
                               ComputeCoarseFissionProduction(coarse_phi_old_),
                               coarse_phi_new_,
                               coarse_production_old);
  const double coarse_phi_tolerance =
    std::min(coarse_pi_phi_tol, coarse_pi_relative_tol * transport_coarse_change);
  const double coarse_k_tolerance =
    std::min(coarse_pi_k_tol, coarse_pi_relative_tol * transport_coarse_change);
  int coarse_pi_iterations = 0;
  bool coarse_pi_converged = false;
  ConvergenceRateEstimator rate_estimator;
  int coarse_ksp_iterations = 0;
  bool coarse_linear_solves_converged = true;
  t0 = NowSeconds();
  // Wielandt-shifted power iteration: each iterate solves (A - F/k_s) x = (1/k - 1/k_s) F phi,
  // which converges by (1/k_0 - 1/k_s) / (1/k_1 - 1/k_s) per iteration instead of the coarse
  // dominance ratio k_1 / k_0. It converges to the fundamental mode only if k_s is above the
  // fundamental eigenvalue k_0; with k_s below k_0 it converges to another mode. So the solve
  // starts unshifted (k_s infinite), which converges to the fundamental mode from any start, and
  // takes the shift k_s = k (1 + wielandt_shift) once the extrapolated error of the k estimate,
  // its change / (1 - rho), is below half the margin, which puts k_s above k_0. The estimate of a
  // poor transport iterate, as in early updates, can be far from k_0. The shift is then raised if
  // the estimate comes within half the margin, so that the shifted operator stays nonsingular.
  double k_shift = std::numeric_limits<double>::infinity();
  int shifted_iterations = 0;
  if (not coarse_eigensolve_failed)
    SetShiftedOperator(k_shift);
  for (int m = 0; m < coarse_pi_max_its and not coarse_eigensolve_failed; ++m)
  {
    const double k_source = 1.0 / (1.0 / k_eff - 1.0 / k_shift);
    // Held values scale with the production of the iterate relative to the transport iterate,
    // which keeps the coarse problem homogeneous in the iterate.
    coarse_phi =
      SolveCoarseSystem(k_source, coarse_phi_m, coarse_production_old / initial_coarse_production);
    coarse_ksp_iterations += last_ksp_iterations_;
    coarse_linear_solves_converged = coarse_linear_solves_converged and last_ksp_converged_;
    coarse_pi_iterations = m + 1;
    // The correction is skipped once any coarse linear solve fails, so further coarse power
    // iterations would be wasted. The converged reason is the same on every rank.
    if (not last_ksp_converged_)
      break;
    const double coarse_production_new = ComputeCoarseFissionProduction(coarse_phi);
    if (not std::isfinite(coarse_production_new) or coarse_production_new <= 0.0)
    {
      coarse_eigensolve_failed = true;
      break;
    }
    const double k_eff_new =
      1.0 / (1.0 / k_shift +
             (1.0 / k_eff - 1.0 / k_shift) * coarse_production_old / coarse_production_new);
    if (not std::isfinite(k_eff_new) or k_eff_new <= 0.0)
    {
      coarse_eigensolve_failed = true;
      break;
    }
    const double k_change = std::fabs(k_eff_new / k_eff - 1.0);
    // The eigenvalue settles much faster than the flux shape, and at a high dominance ratio rho
    // one power iteration changes the shape by only about (1 - rho) times its error. Stopping on
    // either change alone ends the coarse solve while the fundamental-mode shape is still
    // unconverged, so the change is compared with the tolerance scaled by (1 - rho), with rho
    // estimated from successive changes; the first change has no rate estimate yet.
    const double rate = rate_estimator.Update(NormalizedCoarseFluxChange(
      coarse_phi_m, coarse_production_old, coarse_phi, coarse_production_new));

    k_eff = k_eff_new;
    coarse_phi_m = coarse_phi;
    coarse_production_old = coarse_production_new;
    if (std::isfinite(k_shift))
      ++shifted_iterations;

    // Unshifted iterations count toward convergence too; shifted ones need one after the switch.
    if ((std::isinf(k_shift) ? m > 0 : shifted_iterations > 1) and k_change < coarse_k_tolerance and
        rate_estimator.Converged(coarse_phi_tolerance))
    {
      coarse_pi_converged = true;
      break;
    }
    if (std::isinf(k_shift))
    {
      // The rate estimate needs three iterates.
      if (m >= 2 and k_change < 0.5 * wielandt_shift * (1.0 - rate))
      {
        k_shift = k_eff * (1.0 + wielandt_shift);
        SetShiftedOperator(k_shift);
        // The shifted iteration converges at another rate.
        rate_estimator = ConvergenceRateEstimator();
      }
    }
    else if (k_eff > k_shift / (1.0 + 0.5 * wielandt_shift))
    {
      k_shift = k_eff * (1.0 + wielandt_shift);
      SetShiftedOperator(k_shift);
    }
  }
  // The correction and the k handed to the next sweep depend on the amplitude of the coarse
  // solution. The unshifted power iteration ends with production P_0 k / k_0, so that the
  // corrected flux reproduces the coarse k. The shifted iteration scales differently on the way,
  // so restore that amplitude.
  if (coarse_linear_solves_converged and not coarse_eigensolve_failed)
  {
    const double target = initial_coarse_production * k_eff / transport_k_eff;
    const double scale = target / coarse_production_old;
    for (auto& value : coarse_phi)
      value *= scale;
  }
  // A coarse eigensolve that has not converged gives no valid correction. Its iterate need not be
  // near the fundamental mode: in adjoint mode the coarse operator lacks the diagonal dominance in
  // energy that makes it an M-matrix whatever the closure, and the closure of a poor early
  // transport iterate can give it no positive fundamental mode at all.
  if (not coarse_pi_converged and coarse_linear_solves_converged and not coarse_eigensolve_failed)
  {
    coarse_eigensolve_failed = true;
    if (not warned_coarse_pi_limit_)
    {
      warned_coarse_pi_limit_ = true;
      log.Log0Warning() << no_wrap << "CMFD coarse power iteration stopped after "
                        << coarse_pi_max_its << " iterations before meeting its tolerances (flux "
                        << coarse_phi_tolerance << ", k " << coarse_k_tolerance << ") at update "
                        << outer_iteration_
                        << ", so its correction is skipped. This warning is not repeated.";
    }
  }
  const double coarse_pi_time = NowSeconds() - t0;
  const double coarse_phi_change = RelativeL2Difference(coarse_phi, coarse_phi_new_);
  const double coarse_k_eff = k_eff;

  t0 = NowSeconds();
  FluxUpdateDiagnostics correction_diagnostics;
  bool skipped_correction = false;
  std::string correction_skip_reason;
  if (not coarse_linear_solves_converged)
  {
    skipped_correction = true;
    correction_skip_reason = "coarse_linear_solve_not_converged";
    correction_diagnostics = AnalyzeFluxUpdate(phi_new_local_, transport_k_eff);
    if (verbose_ or outer_iteration_ < 5 or outer_iteration_ % 25 == 0)
      log.Log0Warning() << no_wrap
                        << "CMFD correction skipped because a coarse linear solve did not "
                           "converge.";
  }
  else if (coarse_eigensolve_failed)
  {
    skipped_correction = true;
    correction_skip_reason = "coarse_eigensolve_failed";
    correction_diagnostics = AnalyzeFluxUpdate(phi_new_local_, transport_k_eff);
    if (verbose_ or outer_iteration_ < 5 or outer_iteration_ % 25 == 0)
      log.Log0Warning() << no_wrap
                        << "CMFD correction skipped because the coarse eigensolve did not "
                           "converge to a positive finite production and k.";
  }
  else
    k_eff = ApplyFluxCorrection(phi_new_local_,
                                coarse_phi,
                                k_eff / transport_k_eff,
                                cmfd_reference_fission_production,
                                transport_k_eff,
                                correction_diagnostics,
                                skipped_correction,
                                correction_skip_reason);
  if (skipped_correction)
    k_eff = raw_transport_k_eff;
  last_correction_skipped_ = skipped_correction;
  ++(skipped_correction ? num_skipped_corrections_ : num_applied_corrections_);
  consecutive_skipped_corrections_ = skipped_correction ? consecutive_skipped_corrections_ + 1 : 0;
  if (consecutive_skipped_corrections_ >= persistent_skip_warning_count and
      not warned_persistent_skips_)
  {
    warned_persistent_skips_ = true;
    std::ostringstream warning;
    warning << "CMFD skipped the last " << consecutive_skipped_corrections_
            << " corrections (last reason: " << correction_skip_reason
            << "). Power iteration is proceeding without CMFD acceleration.";
    if (correction_skip_reason == "coarse_linear_solve_not_converged")
      warning << " Configure a stronger coarse linear solver with petsc_options, or reduce the "
                 "coarse system with a coarser coarse_mesh or fewer coarse_groups.";
    else
      warning << " Use smaller coarse cells.";
    log.Log0Warning() << no_wrap << warning.str();
  }
  last_correction_skip_reason_ = skipped_correction ? correction_skip_reason : std::string();
  const double accelerated_k_change = std::fabs(k_eff - solver_->GetEigenvalue()) / k_eff;
  // An update whose correction was skipped is unaccelerated, so its changes do not measure the
  // convergence of the accelerated iteration. The owning solver's k and fission-source tests
  // decide convergence otherwise.
  last_update_allows_convergence_ = not skipped_correction;
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
    LogMetric("skipped", skipped_correction ? 1 : 0);
    LogMetric("min_scalar_flux", correction_diagnostics.min_scalar_flux);
    LogMetric("nonfinite_flux", correction_diagnostics.has_nonfinite_flux ? 1 : 0);
    LogMetric("invalid_k", correction_diagnostics.has_invalid_k ? 1 : 0);
    LogMetric("transport_k_eff", transport_k_eff);
    LogMetric("coarse_k_eff", coarse_k_eff);
    LogMetric("applied_k_eff", k_eff);
    LogMetric("coarse_phi_change", coarse_phi_change);
    LogMetric("closure_fallbacks", num_closure_fallbacks_);
    LogMetric("kept_coarse_values", last_kept_coarse_values_);
  }

  if (verbose_ and (outer_iteration_ < 5 or outer_iteration_ % 25 == 0))
  {
    const auto LogMetric = [this](const std::string& metric, const auto value)
    { LogCMFDAcceleratorMetric("timing", outer_iteration_, metric, MetricValue(value)); };
    LogMetric("transport_time", MaxAcrossRanks(transport_time));
    LogMetric("restrict_old_time", MaxAcrossRanks(last_restrict_old_time_));
    LogMetric("restrict_new_time", MaxAcrossRanks(restrict_new_time));
    LogMetric("homogenization_time", MaxAcrossRanks(homogenization_time));
    LogMetric("coarse_flux_cache_time", MaxAcrossRanks(coarse_flux_cache_time));
    LogMetric("face_current_cache_time", MaxAcrossRanks(face_current_cache_time));
    LogMetric("assemble_time", MaxAcrossRanks(assemble_time));
    LogMetric("residual_time", MaxAcrossRanks(residual_time));
    LogMetric("fission_time", MaxAcrossRanks(fission_time));
    LogMetric("coarse_pi_time", MaxAcrossRanks(coarse_pi_time));
    LogMetric("update_time", MaxAcrossRanks(update_time));
    LogMetric("post_time", MaxAcrossRanks(post_time));
    LogMetric("coarse_pi_iterations", coarse_pi_iterations);
    LogMetric("coarse_pi_converged", coarse_pi_converged ? 1 : 0);
    LogMetric("coarse_ksp_iterations", coarse_ksp_iterations);
    LogMetric("coarse_linear_solves_converged", coarse_linear_solves_converged ? 1 : 0);
    LogMetric("operator_balance_residual", transport_balance_residual.relative_l2);
    LogMetric("transport_current_balance_residual", transport_current_balance_residual.relative_l2);
    LogMetric("accelerated_k_change", accelerated_k_change);
  }

  ++outer_iteration_;
  return k_eff;
}

std::string
CMFDAcceleration::GetPowerIterationConvergenceInfo() const
{
  std::ostringstream out;
  if (last_balance_residual_known_)
    out << "CMFD transport-current balance residual = " << std::scientific << std::setprecision(6)
        << last_transport_current_balance_relative_l2_ << std::defaultfloat;
  if (last_correction_skipped_)
  {
    if (last_balance_residual_known_)
      out << ", ";
    out << "CMFD correction = skipped";
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

void
CMFDAcceleration::InitializeHomogenization()
{
  const auto& grid = *do_problem_.GetGrid();
  const auto& unit_cell_matrices = do_problem_.GetUnitCellMatrices();
  const auto& xs_map = do_problem_.GetBlockID2XSMap();
  const auto num_groups = static_cast<std::size_t>(num_groups_);

  // Each local fine cell contributes to one (coarse-cell owner, coarse cell, block) record. The
  // record pattern depends only on the meshes, so it is exchanged once here and every later
  // homogenization sends only flux integrals.
  using RecordKey = std::tuple<int, uint64_t, unsigned int>;
  std::map<RecordKey, double> record_volumes;
  std::vector<RecordKey> fine_cell_keys(grid.GetLocalCellCount());
  std::vector<bool> has_record(grid.GetLocalCellCount(), false);
  std::size_t local_invalid_volumes = 0;
  for (const auto& membership : coarse_mesh_.LocalFineCellMemberships())
  {
    const auto& cell = grid.GetGlobalCell(membership.fine_cell_id);
    const auto& intV_shapeI = unit_cell_matrices[cell.local_id].intV_shapeI;
    const double volume = std::accumulate(intV_shapeI.begin(), intV_shapeI.end(), 0.0);
    if (not std::isfinite(volume) or volume <= 0.0)
      ++local_invalid_volumes;
    const RecordKey key{
      membership.coarse_cell_partition_id, membership.coarse_cell_id, cell.block_id};
    fine_cell_keys[cell.local_id] = key;
    has_record[cell.local_id] = true;
    record_volumes[key] += volume;
  }

  std::map<RecordKey, std::size_t> record_index;
  homogenization_send_counts_.clear();
  std::map<int, std::vector<uint64_t>> pid_keys;
  std::map<int, std::vector<double>> pid_volumes;
  for (const auto& [key, volume] : record_volumes)
  {
    const auto pid = std::get<0>(key);
    record_index[key] = homogenization_send_counts_[pid]++;
    pid_keys[pid].push_back(std::get<1>(key));
    pid_keys[pid].push_back(std::get<2>(key));
    pid_volumes[pid].push_back(volume);
  }
  homogenization_send_records_.assign(grid.GetLocalCellCount(), {0, 0});
  std::size_t local_unassigned_cells = 0;
  for (std::size_t local_id = 0; local_id < fine_cell_keys.size(); ++local_id)
  {
    if (not has_record[local_id])
    {
      ++local_unassigned_cells;
      continue;
    }
    const auto& key = fine_cell_keys[local_id];
    homogenization_send_records_[local_id] = {std::get<0>(key), record_index.at(key)};
  }

  // Each iteration sends num_groups values per record (and the key exchange below sends two).
  // MapAllToAll sums the counts it receives in an int, so the record counts are exchanged first
  // and both the send and the receive sizes are checked before any keys are exchanged.
  const auto values_per_record = std::max<std::size_t>(num_groups, 2);
  std::vector<std::size_t> send_record_counts(opensn::mpi_comm.size(), 0);
  std::size_t total_send_records = 0;
  for (const auto& [pid, count] : homogenization_send_counts_)
  {
    send_record_counts[pid] = count;
    total_send_records += count;
  }
  std::vector<std::size_t> receive_record_counts(opensn::mpi_comm.size(), 0);
  opensn::mpi_comm.all_to_all(send_record_counts, receive_record_counts);
  std::size_t total_receive_records = 0;
  for (const auto count : receive_record_counts)
    total_receive_records += count;
  ThrowIfAnyRankCounted(
    {{local_invalid_volumes, "fine cells have a non-positive or non-finite volume"},
     {local_unassigned_cells, "local fine cells belong to no coarse cell"},
     {FitsMPICount(total_send_records * values_per_record) ? 0 : 1,
      "ranks would send more homogenization values than an MPI count can hold; use more ranks"},
     {FitsMPICount(total_receive_records * values_per_record) ? 0 : 1,
      "ranks would receive more homogenization values than an MPI count can hold; use more "
      "ranks"}});

  const auto received_keys = MapAllToAll(pid_keys);
  const auto received_volumes = MapAllToAll(pid_volumes);

  homogenized_xs_.assign(coarse_mesh_.NumLocalCells(), {});
  homogenization_receive_targets_.clear();
  std::size_t local_invalid_records = 0;
  for (const auto& [pid, keys] : received_keys)
  {
    const auto volumes_it = received_volumes.find(pid);
    if (volumes_it == received_volumes.end() or keys.size() != 2 * volumes_it->second.size())
    {
      ++local_invalid_records;
      continue;
    }
    auto& targets = homogenization_receive_targets_[pid];
    for (std::size_t r = 0; r < volumes_it->second.size(); ++r)
    {
      const auto coarse_cell_id = keys[2 * r];
      const auto block_id = static_cast<unsigned int>(keys[2 * r + 1]);
      if (not coarse_mesh_.HasLocalCoarseCell(coarse_cell_id) or
          xs_map.find(block_id) == xs_map.end())
      {
        ++local_invalid_records;
        continue;
      }
      const auto local_id = coarse_mesh_.LocalCellFromGlobalID(coarse_cell_id).local_id;
      auto& materials = homogenized_xs_[local_id].materials;
      auto material_it = std::find_if(materials.begin(),
                                      materials.end(),
                                      [block_id](const CoarseCellMaterial& material)
                                      { return material.block_id == block_id; });
      if (material_it == materials.end())
      {
        materials.push_back({block_id, 0.0, std::vector<double>(num_groups, 0.0)});
        material_it = std::prev(materials.end());
      }
      material_it->volume += volumes_it->second[r];
      targets.emplace_back(local_id,
                           static_cast<std::size_t>(std::distance(materials.begin(), material_it)));
    }
  }

  std::size_t local_empty_coarse_cells = 0;
  for (auto& coarse_xs : homogenized_xs_)
  {
    local_empty_coarse_cells += coarse_xs.materials.empty() ? 1 : 0;
    coarse_xs.removal.assign(num_coarse_groups_, 0.0);
    coarse_xs.diffusion.assign(num_coarse_groups_, 0.0);
  }
  ThrowIfAnyRankCounted(
    {{local_invalid_records, "received records reference an unknown coarse cell or block"},
     {local_empty_coarse_cells, "coarse cells received no fine-cell records"}});
}

void
CMFDAcceleration::BuildHomogenizedXS(const std::vector<double>& phi, const bool use_flux)
{
  // use_flux, not the size of phi, selects the weights: a rank without cells has an empty phi but
  // must still join the exchange and the reduction below.
  OpenSnInvalidArgumentIf(use_flux and phi.size() != do_problem_.GetPhiNewLocal().size(),
                          "CMFD homogenization scalar flux vector size mismatch.");

  const auto num_groups = static_cast<std::size_t>(num_groups_);
  // The getters condense the scattering and production matrices on demand. They use the cross
  // sections and mode of this assembly, so that all four describe one operator even after the
  // problem's cross sections are replaced.
  assembled_xs_map_ = do_problem_.GetBlockID2XSMap();
  assembled_adjoint_ = do_problem_.GetOptions().adjoint;
  assembled_delayed_production_ = UseDelayedNeutronProduction(do_problem_);
  const bool adjoint = assembled_adjoint_;
  const auto& xs_map = assembled_xs_map_;
  for (auto& coarse_xs : homogenized_xs_)
    for (auto& material : coarse_xs.materials)
      std::fill(material.weights.begin(), material.weights.end(), 0.0);

  // Integrate the scalar flux of each (coarse cell, block) region and gather the integrals on the
  // coarse-cell owner. Without a flux every coarse group uses volume weights below.
  std::size_t local_invalid_buffers = 0;
  if (use_flux)
  {
    const auto& grid = *do_problem_.GetGrid();
    const auto& transport_views = do_problem_.GetCellTransportViews();
    const auto& unit_cell_matrices = do_problem_.GetUnitCellMatrices();
    std::map<int, std::vector<double>> pid_values;
    for (const auto& [pid, count] : homogenization_send_counts_)
      pid_values[pid].assign(count * num_groups, 0.0);
    for (const auto& cell : grid.GetLocalCells())
    {
      const auto& [pid, record] = homogenization_send_records_[cell->local_id];
      auto& integrals = pid_values.at(pid);
      const auto offset = record * num_groups;
      const auto& view = transport_views[cell->local_id];
      const auto& intV_shapeI = unit_cell_matrices[cell->local_id].intV_shapeI;
      for (int i = 0; i < view.GetNumNodes(); ++i)
      {
        const auto phi_map = view.MapDOF(i, 0, first_group_);
        for (std::size_t g = 0; g < num_groups; ++g)
          integrals[offset + g] += phi[phi_map + g] * intV_shapeI(i);
      }
    }

    const auto received_values = MapAllToAll(pid_values);
    for (const auto& [pid, values] : received_values)
    {
      const auto targets_it = homogenization_receive_targets_.find(pid);
      if (targets_it == homogenization_receive_targets_.end() or
          values.size() != targets_it->second.size() * num_groups)
      {
        ++local_invalid_buffers;
        continue;
      }
      for (std::size_t r = 0; r < targets_it->second.size(); ++r)
      {
        const auto& [local_id, material_index] = targets_it->second[r];
        auto& weights = homogenized_xs_[local_id].materials[material_index].weights;
        for (std::size_t g = 0; g < num_groups; ++g)
          weights[g] += values[r * num_groups + g];
      }
    }
  }

  // Normalize each coarse group's flux integrals into weights. Signed weights preserve the
  // transport reaction rates exactly, so the CMFD fixed point matches the transport solution;
  // volume weights are used only where the signed coarse-group flux is zero, non-finite, or lost
  // to cancellation.
  std::size_t local_fallback_groups = 0;
  for (auto& coarse_xs : homogenized_xs_)
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      const auto g_begin = FineGroupBegin(cg);
      const auto g_end = FineGroupEnd(cg);
      double signed_sum = 0.0;
      double absolute_sum = 0.0;
      for (const auto& material : coarse_xs.materials)
        for (auto g = g_begin; g < g_end; ++g)
        {
          signed_sum += material.weights[g];
          absolute_sum += std::fabs(material.weights[g]);
        }
      const bool flux_weighted =
        use_flux and std::isfinite(absolute_sum) and
        std::fabs(signed_sum) > homogenization_cancellation_tolerance * absolute_sum;
      if (flux_weighted)
        for (auto& material : coarse_xs.materials)
          for (auto g = g_begin; g < g_end; ++g)
            material.weights[g] /= signed_sum;
      else
      {
        local_fallback_groups += use_flux ? 1 : 0;
        double volume_sum = 0.0;
        for (const auto& material : coarse_xs.materials)
          volume_sum += material.volume * static_cast<double>(g_end - g_begin);
        for (auto& material : coarse_xs.materials)
          for (auto g = g_begin; g < g_end; ++g)
            material.weights[g] = material.volume / volume_sum;
      }

      // Removal is total interaction minus scattering that stays within the coarse group.
      double removal = 0.0;
      for (const auto& material : coarse_xs.materials)
      {
        const auto& xs = *xs_map.at(material.block_id);
        const auto& sigma_t = xs.GetSigmaTotal();
        const auto& scattering = xs.GetTransferMatrix(0, adjoint);
        for (auto g = g_begin; g < g_end; ++g)
        {
          removal += sigma_t[g] * material.weights[g];
          for (const auto& [_, gp, sigma_s] : scattering.Row(g))
            if (gp >= g_begin and gp < g_end)
              removal -= sigma_s * material.weights[gp];
        }
      }
      coarse_xs.removal[cg] = removal;

      // The diffusion coefficient only sets the coupling that the current correction completes,
      // so it uses weight magnitudes, which keep it positive. Each fine group's transport cross
      // section is homogenized in space before it is inverted.
      double diffusion_numerator = 0.0;
      double diffusion_denominator = 0.0;
      for (auto g = g_begin; g < g_end; ++g)
      {
        double group_weight = 0.0;
        double transport_rate = 0.0;
        for (const auto& material : coarse_xs.materials)
        {
          const double weight = std::fabs(material.weights[g]);
          group_weight += weight;
          transport_rate += xs_map.at(material.block_id)->GetSigmaTransport(adjoint)[g] * weight;
        }
        if (group_weight == 0.0)
          continue;
        const double sigma_tr = transport_rate / group_weight;
        const double diffusion = sigma_tr < MultiGroupXS::TRANSPORT_XS_TOLERANCE
                                   ? MultiGroupXS::MAX_DIFFUSION_COEFFICIENT
                                   : 1.0 / (3.0 * sigma_tr);
        diffusion_numerator += group_weight * diffusion;
        diffusion_denominator += group_weight;
      }
      coarse_xs.diffusion[cg] = diffusion_numerator / diffusion_denominator;
    }

  if (not use_flux)
    return;

  std::vector<std::size_t> local_counts = {local_fallback_groups, local_invalid_buffers};
  std::vector<std::size_t> counts(local_counts.size(), 0);
  opensn::mpi_comm.all_reduce(local_counts, counts, mpi::op::sum<std::size_t>());
  if (counts[1] > 0)
    throw std::logic_error("CMFD homogenization: " + std::to_string(counts[1]) +
                           " received flux-integral buffers do not match the exchange pattern.");
  if (counts[0] > 0 and not warned_homogenization_fallback_)
  {
    log.Log0Warning() << "CMFD used volume weights for " << counts[0]
                      << " coarse-cell energy groups whose scalar-flux sum was zero, non-finite, "
                         "or cancellation dominated. Their coarse rates do not match the "
                         "transport rates.";
    warned_homogenization_fallback_ = true;
  }
}

void
CMFDAcceleration::CheckAssembledGroupCount(const MultiGroupXS& xs,
                                           const unsigned int block_id) const
{
  // The cross-section objects of the last assembly are shared with the problem, so a caller can
  // change them in place after a solve. The condensation loops over the solver's groups.
  OpenSnLogicalErrorIf(xs.GetNumGroups() != num_groups_,
                       GetName() + ": the cross sections of block " + std::to_string(block_id) +
                         " no longer have the group count of the last CMFD assembly.");
}

void
CMFDAcceleration::ComputeCoarseScattering(const CoarseCellXS& coarse_xs,
                                          std::vector<double>& matrix) const
{
  const auto n = static_cast<std::size_t>(num_coarse_groups_);
  const bool adjoint = assembled_adjoint_;
  const auto& xs_map = assembled_xs_map_;
  matrix.assign(n * n, 0.0);
  for (const auto& material : coarse_xs.materials)
  {
    CheckAssembledGroupCount(*xs_map.at(material.block_id), material.block_id);
    const auto& scattering = xs_map.at(material.block_id)->GetTransferMatrix(0, adjoint);
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      const auto row = static_cast<std::size_t>(g / group_aggregation_size_) * n;
      for (const auto& [_, gp, sigma_s] : scattering.Row(g))
        matrix[row + gp / group_aggregation_size_] += sigma_s * material.weights[gp];
    }
  }
}

bool
CMFDAcceleration::ComputeCoarseProduction(const CoarseCellXS& coarse_xs,
                                          std::vector<double>& matrix) const
{
  const auto n = static_cast<std::size_t>(num_coarse_groups_);
  const bool adjoint = assembled_adjoint_;
  const auto& xs_map = assembled_xs_map_;
  matrix.assign(n * n, 0.0);
  bool fissionable = false;
  for (const auto& material : coarse_xs.materials)
  {
    const auto& xs = *xs_map.at(material.block_id);
    CheckAssembledGroupCount(xs, material.block_id);
    if (not xs.IsFissionable())
      continue;
    fissionable = true;
    // In adjoint mode the production matrix is transposed, and so is the delayed production.
    const auto& production = xs.GetProductionMatrix(adjoint);
    const bool include_delayed_production = assembled_delayed_production_;
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      const auto row = static_cast<std::size_t>(g / group_aggregation_size_) * n;
      for (unsigned int gp = 0; gp < num_groups_; ++gp)
      {
        double value = production[g][gp];
        if (include_delayed_production)
          value += adjoint ? ComputeDelayedFissionProduction(xs, gp, g)
                           : ComputeDelayedFissionProduction(xs, g, gp);
        matrix[row + gp / group_aggregation_size_] += value * material.weights[gp];
      }
    }
  }
  return fissionable;
}

void
CMFDAcceleration::InitializeLinearSystem()
{
  // Guard against leaking these if Initialize is called multiple times.
  if (ksp_)
    OpenSnPETScCall(KSPDestroy(&ksp_));
  if (rhs_)
    OpenSnPETScCall(VecDestroy(&rhs_));
  if (A_)
    OpenSnPETScCall(MatDestroy(&A_));
  if (A_shift_)
    OpenSnPETScCall(MatDestroy(&A_shift_));

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

  // Accelerations of other problems can have the same name, so the coarse solver keeps its options
  // apart and puts them in the global database only while it uses them.
  petsc_options_db_ = std::make_unique<PETScOptionsDatabase>();
  OpenSnPETScCall(KSPCreate(opensn::mpi_comm, &ksp_));
  OpenSnPETScCall(KSPSetOperators(ksp_, A_, A_));
  OpenSnPETScCall(KSPSetOptionsPrefix(ksp_, GetName().c_str()));
  ConfigureCoarseLinearSolver();
  OpenSnPETScCall(KSPSetTolerances(
    ksp_, coarse_linear_rel_tol, coarse_linear_abs_tol, PETSC_DEFAULT, coarse_linear_max_its));
  // The options are in the global database only while a scope of petsc_options_db_ exists, around
  // the calls of the coarse solver that read them: KSPSetFromOptions, KSPSetUp and KSPSolve.
  coarse_solver_from_options_ = not petsc_options_.empty();
  if (coarse_solver_from_options_)
    petsc_options_db_->Insert(petsc_options_);
  const PETScOptionsDatabase::Scope petsc_options_scope(*petsc_options_db_);
  OpenSnPETScCall(KSPSetFromOptions(ksp_));
}

void
CMFDAcceleration::ConfigureCoarseLinearSolver()
{
  // Estimated memory of an LU factorization with nested-dissection ordering. A mesh-like graph of
  // n coarse cells fills in to about n entries in 1D, 4 n log2(n) in 2D and 6 n^(4/3) in 3D, each
  // a dense G x G block, at about 12 bytes per stored entry for L and U.
  const auto num_cells = static_cast<double>(coarse_mesh_.NumGlobalCells());
  const auto groups = static_cast<double>(num_coarse_groups_);
  const auto dimension = do_problem_.GetGrid()->GetDimension();
  const double cell_fill = dimension <= 1   ? 3.0 * num_cells
                           : dimension == 2 ? 4.0 * num_cells * std::log2(std::max(num_cells, 2.0))
                                            : 6.0 * std::pow(num_cells, 4.0 / 3.0);
  const double lu_bytes = 2.0 * 12.0 * groups * groups * cell_fill;
  const double budget = coarse_direct_bytes_per_rank * static_cast<double>(opensn::mpi_comm.size());
  bool direct = lu_bytes <= budget;

  PC pc = nullptr;
  OpenSnPETScCall(KSPGetPC(ksp_, &pc));
  std::string solver_name;
  if (direct)
  {
    OpenSnPETScCall(KSPSetType(ksp_, KSPPREONLY));
    OpenSnPETScCall(PCSetType(pc, PCLU));
#ifdef PETSC_HAVE_SUPERLU_DIST
    OpenSnPETScCall(PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST));
    solver_name = "LU (SuperLU_DIST)";
#elif defined(PETSC_HAVE_MUMPS)
    OpenSnPETScCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    solver_name = "LU (MUMPS)";
#else
    // PETSc's own LU is sequential.
    direct = opensn::mpi_comm.size() == 1;
    solver_name = "LU (PETSc)";
#endif
  }
  if (not direct)
  {
    // Algebraic multigrid can report convergence on coarse operators that are not M-matrices,
    // which signed transport iterates produce. Additive Schwarz with ILU subdomain solves needs no
    // M-matrix, and the true-residual test in SolveCoarseSystem rejects a failed solve.
    OpenSnPETScCall(KSPSetType(ksp_, KSPFGMRES));
    OpenSnPETScCall(KSPGMRESSetRestart(ksp_, 100));
    OpenSnPETScCall(PCSetType(pc, PCASM));
    OpenSnPETScCall(PCASMSetOverlap(pc, 1));
    solver_name = "FGMRES with additive Schwarz and ILU(0)";
  }
  using_direct_coarse_solver_ = direct;
  log.Log() << no_wrap << program_timer.GetTimeString()
            << " CMFD coarse linear solver: " << solver_name << " (estimated LU memory "
            << lu_bytes / (1024.0 * 1024.0 * 1024.0) << " GiB, budget "
            << coarse_direct_bytes_per_rank / (1024.0 * 1024.0 * 1024.0) << " GiB per rank, "
            << budget / (1024.0 * 1024.0 * 1024.0) << " GiB in total)"
            << (petsc_options_.empty() ? "." : ", then petsc_options.");
}

void
CMFDAcceleration::ConfigureTransportSolve(const unsigned int max_iterations,
                                          const double residual_tolerance)
{
  for (std::size_t gsid = 0; gsid < do_problem_.GetNumGroupsets(); ++gsid)
    do_problem_.SetWGSSolveSettings(gsid, max_iterations, residual_tolerance);
}

double
CMFDAcceleration::CachedCoarseFlux(const std::size_t local_index) const
{
  // A value held in the previous coarse solve stays held, and its neighbors see that from a
  // non-positive cached flux. The closure reads the flux of a free neighbor only.
  const double phi = coarse_phi_new_[local_index];
  if (held_.size() == coarse_phi_new_.size() and held_[local_index])
    return std::min(phi, 0.0);
  return phi;
}

void
CMFDAcceleration::BuildCoarseFluxCache()
{
  coarse_phi_cache_.clear();
  diffusion_coefficient_cache_.clear();
  const bool has_coarse_phi = not coarse_phi_new_.empty();
  const std::size_t values_per_key = has_coarse_phi ? 2 : 1;

  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      const auto key = std::make_tuple(coarse_cell.global_id, cg);
      diffusion_coefficient_cache_[key] = homogenized_xs_[coarse_cell.local_id].diffusion[cg];
      if (has_coarse_phi)
        coarse_phi_cache_[key] = CachedCoarseFlux(offset + cg);
    }
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
    response_values.reserve(values_per_key * requests.size() / 2);

    for (std::size_t r = 0; r < requests.size(); r += 2)
    {
      const auto coarse_cell_gid = requests[r];
      const auto cg = static_cast<unsigned int>(requests[r + 1]);
      OpenSnInvalidArgumentIf(not coarse_mesh_.HasLocalCoarseCell(coarse_cell_gid),
                              "CMFD coarse-flux request references a nonlocal coarse cell.");

      const auto coarse_local_id = coarse_mesh_.LocalCellFromGlobalID(coarse_cell_gid).local_id;
      response_keys.push_back(coarse_cell_gid);
      response_keys.push_back(cg);
      response_values.push_back(homogenized_xs_[coarse_local_id].diffusion[cg]);
      if (has_coarse_phi)
        response_values.push_back(CachedCoarseFlux(coarse_local_id * num_coarse_groups_ + cg));
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
    OpenSnLogicalErrorIf(keys.size() % 2 != 0 or values.size() != values_per_key * keys.size() / 2,
                         "Invalid CMFD coarse-flux response buffer.");

    for (std::size_t r = 0; r < keys.size() / 2; ++r)
    {
      const auto key = std::make_tuple(keys[2 * r], static_cast<unsigned int>(keys[2 * r + 1]));
      diffusion_coefficient_cache_[key] = values[values_per_key * r];
      if (has_coarse_phi)
        coarse_phi_cache_[key] = values[values_per_key * r + 1];
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
      const auto& owner_cell = grid.GetGlobalCell(owner_cell_gid);
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

CMFDAcceleration::OperatorSignCounts
CMFDAcceleration::AssembleOperator(const bool from_transport_iterate)
{
  OpenSnPETScCall(MatZeroEntries(A_));

  // A coarse flux this small relative to the largest flux of its coarse group is zero for the
  // closure, whose current corrections divide by it. The flag, not the local coarse flux, decides
  // the collectives: a rank may own no coarse cells.
  std::vector<double> zero_flux(num_coarse_groups_, 0.0);
  if (from_transport_iterate)
  {
    std::vector<double> local_max(num_coarse_groups_, 0.0);
    for (std::size_t i = 0; i < coarse_phi_new_.size(); ++i)
      local_max[i % num_coarse_groups_] =
        std::max(local_max[i % num_coarse_groups_], std::fabs(coarse_phi_new_[i]));
    opensn::mpi_comm.all_reduce(local_max.data(),
                                static_cast<int>(num_coarse_groups_),
                                zero_flux.data(),
                                mpi::op::max<double>());
    for (auto& value : zero_flux)
      value *= zero_flux_fraction;
  }

  // A coarse value whose restricted transport flux is not positive is held at that flux times the
  // amplitude of the coarse iterate (see AssembleRHS). The partial-current closure divides by the
  // coarse flux, so with a negative flux it gives the operator entries of the wrong sign, and the
  // coarse solution near such values is unreliable. A free value next to a held one models its
  // outgoing partial current by the ratio to its own flux and receives the held value's outgoing
  // partial current as a source, scaled like the held value. At the transport solution both
  // reproduce the transport partial currents, so the fixed point is unchanged.
  // A value is held only while its restricted transport flux is not positive. A held value fixes
  // a source in the coarse eigenproblem, so a value held after its flux has become positive again
  // keeps an early iterate's flux in every later coarse solve.
  held_.assign(LocalUnknownCount(), 0);
  held_source_.assign(LocalUnknownCount(), 0.0);
  if (from_transport_iterate)
    for (std::size_t i = 0; i < coarse_phi_new_.size(); ++i)
      if (coarse_phi_new_[i] <= zero_flux[i % num_coarse_groups_])
      {
        held_[i] = 1;
        held_source_[i] = coarse_phi_new_[i];
      }

  OperatorSignCounts sign_counts;
  std::size_t local_fallbacks = 0;
  std::vector<double> scattering;
  // Face couplings to one neighbor can arrive from several coarse faces, so the row's
  // off-diagonal entries are summed here before their signs are counted.
  std::vector<std::pair<int64_t, double>> row_face_entries;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = homogenized_xs_.at(coarse_cell.local_id);
    ComputeCoarseScattering(xs, scattering);

    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      const auto row = MapDOF(coarse_cell.global_id, cg);
      const auto local_index = coarse_cell.local_id * num_coarse_groups_ + cg;
      if (held_[local_index])
      {
        OpenSnPETScCall(MatSetValue(A_, row, row, 1.0, ADD_VALUES));
        continue;
      }
      double diag = xs.removal[cg] * coarse_cell.volume;
      row_face_entries.clear();

      for (std::size_t f = 0; f < coarse_cell.faces.size(); ++f)
      {
        const auto& face = coarse_cell.faces[f];
        if (face.has_neighbor)
        {
          // Diffusion coupling (D-tilde). Before a transport iterate exists, or where a coarse
          // flux is zero, the face is closed by this coupling alone.
          const double coeff = InterfaceDiffusionCoefficient(
            HalfTransmissibility(xs.diffusion[cg], coarse_cell.centroid, face.centroid, face.area),
            HalfTransmissibility(diffusion_coefficient_cache_.at({face.neighbor_id, cg}),
                                 face.neighbor_centroid,
                                 face.centroid,
                                 face.area));
          double face_diag = coeff;
          double face_offdiag = -coeff;
          if (from_transport_iterate and
              coarse_phi_cache_.at(std::make_tuple(face.neighbor_id, cg)) <= zero_flux[cg])
          {
            // The neighbor is held: the outgoing partial current is the ratio form, and the
            // incoming one is a source.
            const auto phi_owner = coarse_phi_new_[local_index];
            const auto [owner_partial_current, neighbor_partial_current] =
              ComputePartialOutwardCurrents(coarse_cell, f, cg);
            diag += owner_partial_current / phi_owner;
            held_source_[local_index] += neighbor_partial_current;
            continue;
          }
          if (from_transport_iterate)
          {
            const auto phi_owner = coarse_phi_new_[local_index];
            const auto phi_neighbor = coarse_phi_cache_.at(std::make_tuple(face.neighbor_id, cg));
            // pCMFD closure: each partial current carries half of the diffusion coupling,
            // J+ = D/2 (phi_owner - phi_neighbor) + Dhat+ phi_owner and
            // J- = D/2 (phi_neighbor - phi_owner) + Dhat- phi_neighbor, with Dhat+- chosen so
            // that J+- reproduce the transport partial currents.
            bool closed = false;
            if (std::fabs(phi_owner) > zero_flux[cg] and std::fabs(phi_neighbor) > zero_flux[cg])
            {
              const auto [owner_partial_current, neighbor_partial_current] =
                ComputePartialOutwardCurrents(coarse_cell, f, cg);
              const double half_coeff = 0.5 * coeff;
              const double owner_dhat =
                (owner_partial_current - half_coeff * (phi_owner - phi_neighbor)) / phi_owner;
              const double neighbor_dhat =
                (neighbor_partial_current - half_coeff * (phi_neighbor - phi_owner)) / phi_neighbor;
              if (std::isfinite(owner_dhat) and std::isfinite(neighbor_dhat))
              {
                face_diag = coeff + owner_dhat;
                face_offdiag = -coeff - neighbor_dhat;
                closed = true;
              }
            }
            if (not closed)
              ++local_fallbacks;
          }

          diag += face_diag;
          const auto col = MapDOF(face.neighbor_id, cg);
          OpenSnPETScCall(MatSetValue(A_, row, col, face_offdiag, ADD_VALUES));
          const auto entry =
            std::find_if(row_face_entries.begin(),
                         row_face_entries.end(),
                         [col](const std::pair<int64_t, double>& e) { return e.first == col; });
          if (entry == row_face_entries.end())
            row_face_entries.emplace_back(col, face_offdiag);
          else
            entry->second += face_offdiag;
        }
        else
        {
          const auto bc_it = do_problem_.GetBoundaryDefinitions().find(face.neighbor_id);
          const bool reflective = bc_it != do_problem_.GetBoundaryDefinitions().end() and
                                  bc_it->second.type == LBSBoundaryType::REFLECTING;
          if (not reflective)
          {
            // Marshak coupling until a transport iterate exists. Afterwards the coefficient
            // reproduces the transport outgoing current.
            double boundary_coeff = 0.5 * face.area;
            if (from_transport_iterate)
            {
              const auto phi_owner =
                coarse_phi_new_[coarse_cell.local_id * num_coarse_groups_ + cg];
              if (std::fabs(phi_owner) > zero_flux[cg])
                boundary_coeff = ComputeOutwardCurrent(coarse_cell, f, cg) / phi_owner;
              else
                ++local_fallbacks;
            }
            diag += boundary_coeff;
          }
        }
      }

      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        if (cgp != cg)
        {
          const auto sigma_s = scattering[cg * num_coarse_groups_ + cgp];
          const auto col = MapDOF(coarse_cell.global_id, cgp);
          OpenSnPETScCall(MatSetValue(A_, row, col, -sigma_s * coarse_cell.volume, ADD_VALUES));
          if (-sigma_s * coarse_cell.volume > 0.0)
            ++sign_counts.positive_off_diagonals;
        }

      OpenSnPETScCall(MatSetValue(A_, row, row, diag, ADD_VALUES));
      if (not(diag > 0.0))
        ++sign_counts.nonpositive_diagonals;
      for (const auto& entry : row_face_entries)
        if (entry.second > 0.0)
          ++sign_counts.positive_off_diagonals;
    }
  }

  OpenSnPETScCall(MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY));
  OpenSnPETScCall(MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY));
  OpenSnPETScCall(KSPSetOperators(ksp_, A_, A_));

  if (from_transport_iterate)
  {
    std::size_t global_fallbacks = 0;
    opensn::mpi_comm.all_reduce(local_fallbacks, global_fallbacks, mpi::op::sum<std::size_t>());
    if (global_fallbacks > 0 and num_closure_fallbacks_ == 0)
      log.Log0Warning() << no_wrap << "CMFD closed " << global_fallbacks
                        << " coarse faces by the diffusion coupling alone at update "
                        << outer_iteration_
                        << " because a coarse flux next to them is zero. The coarse currents "
                           "there do not reproduce the transport currents. This warning is not "
                           "repeated.";
    num_closure_fallbacks_ += global_fallbacks;
  }
  return sign_counts;
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
      const auto& fine_cell = grid.GetGlobalCell(fine_face.cell_id);
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
      const auto& neighbor_cell = grid.GetGlobalCell(*neighbor_id);
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
                                    const std::vector<double>& coarse_source_phi,
                                    const double held_scale) const
{
  OpenSnInvalidArgumentIf(coarse_source_phi.size() !=
                            coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse source scalar flux vector size mismatch.");
  const PETScOptionsDatabase::Scope petsc_options_scope(*petsc_options_db_);

  AssembleRHS(rhs_, k_eff, coarse_source_phi, held_scale);

  Vec x = nullptr;
  OpenSnPETScCall(VecDuplicate(rhs_, &x));
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      OpenSnPETScCall(VecSetValue(x,
                                  MapDOF(coarse_cell.global_id, cg),
                                  coarse_source_phi[coarse_offset + cg],
                                  INSERT_VALUES));
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

  // PETSc reports success for every direct solve and tests the preconditioned residual of
  // iterative ones, so a near-singular coarse operator can pass with a useless solution. The
  // solve is accepted only if the true residual is within true_residual_factor of the requested
  // tolerance. The factor allows for the gap between preconditioned and true residual norms.
  Mat ksp_operator = nullptr;
  OpenSnPETScCall(KSPGetOperators(ksp_, &ksp_operator, nullptr));
  Vec residual = nullptr;
  OpenSnPETScCall(VecDuplicate(rhs_, &residual));
  OpenSnPETScCall(MatMult(ksp_operator, x, residual));
  OpenSnPETScCall(VecAYPX(residual, -1.0, rhs_));
  PetscReal residual_norm = 0.0;
  PetscReal rhs_norm = 0.0;
  OpenSnPETScCall(VecNorm(residual, NORM_2, &residual_norm));
  OpenSnPETScCall(VecNorm(rhs_, NORM_2, &rhs_norm));
  OpenSnPETScCall(VecDestroy(&residual));
  const double true_residual_tolerance =
    true_residual_factor *
    std::max(coarse_linear_rel_tol * static_cast<double>(rhs_norm), coarse_linear_abs_tol);
  last_ksp_converged_ = ksp_reason > 0 and std::isfinite(static_cast<double>(residual_norm)) and
                        static_cast<double>(residual_norm) <= true_residual_tolerance;

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
                              const std::vector<double>& coarse_source_phi,
                              const double held_scale) const
{
  OpenSnInvalidArgumentIf(coarse_source_phi.size() !=
                            coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse source scalar flux vector size mismatch.");

  OpenSnPETScCall(VecSet(rhs, 0.0));
  // A held value's row is its restricted transport flux, and a free value next to held ones
  // receives their outgoing partial currents, both times the amplitude of the coarse iterate.
  const bool has_held = held_.size() == coarse_source_phi.size();
  if (has_held)
    for (const auto& coarse_cell : coarse_mesh_.LocalCells())
      for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      {
        const auto i = coarse_cell.local_id * num_coarse_groups_ + cg;
        if (held_source_[i] != 0.0)
          OpenSnPETScCall(VecSetValue(
            rhs, MapDOF(coarse_cell.global_id, cg), held_scale * held_source_[i], ADD_VALUES));
      }
  std::vector<double> production;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    if (not ComputeCoarseProduction(homogenized_xs_.at(coarse_cell.local_id), production))
      continue;

    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      if (has_held and held_[coarse_offset + cg])
        continue;
      double rhs_value = 0.0;
      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        rhs_value +=
          production[cg * num_coarse_groups_ + cgp] * coarse_source_phi[coarse_offset + cgp];
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

  AssembleRHS(rhs, k_eff, coarse_source_phi, 1.0);
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
          rhs_l2 > 0.0 ? static_cast<double>(residual_l2 / rhs_l2)
                       : std::numeric_limits<double>::infinity()};
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

  double local_max_abs = 0.0;
  std::vector<std::pair<double, double>> residual_rhs;
  residual_rhs.reserve(LocalUnknownCount());

  std::vector<double> scattering;
  std::vector<double> production;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = homogenized_xs_.at(coarse_cell.local_id);
    ComputeCoarseScattering(xs, scattering);
    ComputeCoarseProduction(xs, production);
    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;

    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      const double phi = coarse_phi[coarse_offset + cg];
      double lhs = xs.removal[cg] * coarse_cell.volume * phi;

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
          const auto sigma_s = scattering[cg * num_coarse_groups_ + cgp];
          lhs -= sigma_s * coarse_cell.volume * coarse_phi[coarse_offset + cgp];
        }

      double rhs = 0.0;
      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        rhs += production[cg * num_coarse_groups_ + cgp] * coarse_source_phi[coarse_offset + cgp];
      rhs *= coarse_cell.volume / k_eff;

      const double residual = lhs - rhs;
      local_max_abs = std::max(local_max_abs, std::fabs(residual));
      residual_rhs.emplace_back(residual, rhs);
    }
  }

  double global_max_abs = 0.0;
  mpi_comm.all_reduce(local_max_abs, global_max_abs, mpi::op::max<double>());
  return {global_max_abs, CMFDRelativeBalanceResidual(residual_rhs)};
}

double
CMFDAcceleration::ComputeCoarseFissionProduction(const std::vector<double>& coarse_phi) const
{
  OpenSnInvalidArgumentIf(coarse_phi.size() != coarse_mesh_.NumLocalCells() * num_coarse_groups_,
                          "Coarse scalar flux vector size mismatch.");

  double local_production = 0.0;
  std::vector<double> production;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    if (not ComputeCoarseProduction(homogenized_xs_.at(coarse_cell.local_id), production))
      continue;

    const auto coarse_offset = coarse_cell.local_id * num_coarse_groups_;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        local_production += production[cg * num_coarse_groups_ + cgp] *
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
  for (const auto& cell : grid.GetLocalCells())
  {
    const auto& transport_view = transport_views[cell->local_id];
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
CMFDAcceleration::ApplyFluxCorrection(const std::vector<double>& transport_phi_new,
                                      const std::vector<double>& coarse_phi,
                                      const double coarse_amplitude,
                                      const double old_fission_production,
                                      const double transport_k_eff,
                                      FluxUpdateDiagnostics& diagnostics,
                                      bool& skipped_correction,
                                      std::string& skip_reason)
{
  skipped_correction = false;
  skip_reason.clear();
  const auto transport_diagnostics = AnalyzeFluxUpdate(transport_phi_new, transport_k_eff);
  OpenSnLogicalErrorIf(transport_diagnostics.has_invalid_k,
                       "CMFD input transport update has an invalid k-eigenvalue.");
  OpenSnLogicalErrorIf(transport_diagnostics.has_nonfinite_flux,
                       "CMFD input transport update has non-finite scalar flux values.");

  const auto Skip = [&](const std::string& reason)
  {
    diagnostics = transport_diagnostics;
    phi_new_local_ = transport_phi_new;
    skipped_correction = true;
    skip_reason = reason;
    if (verbose_ or outer_iteration_ < 5 or outer_iteration_ % 25 == 0)
      log.Log0Warning() << no_wrap << "CMFD correction skipped (" << reason
                        << "). Using the unaccelerated transport update.";
    return transport_k_eff;
  };

  // The correction scales the fine flux of each coarse cell and coarse group by a coarse ratio. A
  // negative or non-finite ratio would change the sign of the transport flux it scales, which
  // happens where the transport iterate or the coarse solution has the wrong sign, for example in
  // elongated cells of a discontinuous discretization. Such a coarse value keeps the transport
  // iterate, scaled to the amplitude of the coarse solution, and the other values are corrected.
  // The ratios tend to one at convergence, so the fixed point is unchanged.
  auto admissible_phi = coarse_phi;
  std::size_t local_kept = 0;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      const auto i = coarse_cell.local_id * num_coarse_groups_ + cg;
      if (not CMFDRatioChangesSign(coarse_phi_new_[i], coarse_phi[i]))
        continue;
      if (verbose_ and local_kept < 3)
        log.LogAll() << no_wrap << "CMFD keeps the transport flux at coarse cell "
                     << coarse_cell.global_id << " (" << coarse_cell.centroid.x << ", "
                     << coarse_cell.centroid.y << ", " << coarse_cell.centroid.z
                     << "), coarse group " << cg << ": transport " << coarse_phi_new_[i]
                     << ", coarse " << coarse_phi[i];
      admissible_phi[i] = coarse_phi_new_[i] * coarse_amplitude;
      ++local_kept;
    }
  opensn::mpi_comm.all_reduce(local_kept, last_kept_coarse_values_, mpi::op::sum<std::size_t>());
  if (last_kept_coarse_values_ > 0 and not warned_kept_coarse_values_)
  {
    warned_kept_coarse_values_ = true;
    log.Log0Warning() << no_wrap << "CMFD kept the transport flux in " << last_kept_coarse_values_
                      << " of " << GlobalUnknownCount() << " coarse values at update "
                      << outer_iteration_
                      << " because the restricted transport flux or the coarse solution there is "
                         "negative, so the correction would change or amplify the sign of the "
                         "flux. This happens in early updates and where a discontinuous solution "
                         "dips below zero, for example in reflectors. verbose=True lists where. "
                         "This warning is not repeated.";
  }

  auto candidate_phi_new = transport_phi_new;
  const auto fine_cell_ratios = CMFDComputeFineCellRatios(do_problem_,
                                                          num_groups_,
                                                          group_aggregation_size_,
                                                          coarse_mesh_,
                                                          coarse_phi_new_,
                                                          admissible_phi);
  CMFDApplyFineCellRatios(do_problem_,
                          first_group_,
                          num_groups_,
                          group_aggregation_size_,
                          fine_cell_ratios,
                          false,
                          candidate_phi_new);
  const double candidate_fission_production =
    ComputeFissionProduction(do_problem_, candidate_phi_new);
  const double candidate_k_eff =
    candidate_fission_production / old_fission_production * solver_->GetEigenvalue();
  diagnostics = AnalyzeFluxUpdate(candidate_phi_new, candidate_k_eff);
  if (diagnostics.has_invalid_k)
    return Skip("invalid_k");
  if (diagnostics.has_nonfinite_flux)
    return Skip("nonfinite_flux");
  phi_new_local_ = std::move(candidate_phi_new);
  return candidate_k_eff;
}

void
CMFDAcceleration::SetShiftedOperator(const double k_shift)
{
  OpenSnLogicalErrorIf(not(k_shift > 0.0), "CMFD Wielandt shift must be positive.");
  // An infinite shift is the unshifted operator.
  if (std::isinf(k_shift))
  {
    OpenSnPETScCall(KSPSetOperators(ksp_, A_, A_));
    const PETScOptionsDatabase::Scope petsc_options_scope(*petsc_options_db_);
    OpenSnPETScCall(KSPSetUp(ksp_));
    return;
  }
  if (A_shift_)
    OpenSnPETScCall(MatDestroy(&A_shift_));
  OpenSnPETScCall(MatDuplicate(A_, MAT_COPY_VALUES, &A_shift_));
  OpenSnPETScCall(MatSetOption(A_shift_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

  // The operator couples homogenized coarse cells, so the shift subtracts the homogenized
  // production of each coarse cell, the same matrix the right-hand side uses.
  std::vector<double> production;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    if (not ComputeCoarseProduction(homogenized_xs_.at(coarse_cell.local_id), production))
      continue;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      if (not held_.empty() and held_[coarse_cell.local_id * num_coarse_groups_ + cg])
        continue;
      const auto row = MapDOF(coarse_cell.global_id, cg);
      for (unsigned int cgp = 0; cgp < num_coarse_groups_; ++cgp)
        if (const double value = production[cg * num_coarse_groups_ + cgp]; value != 0.0)
          OpenSnPETScCall(MatSetValue(A_shift_,
                                      row,
                                      MapDOF(coarse_cell.global_id, cgp),
                                      -value * coarse_cell.volume / k_shift,
                                      ADD_VALUES));
    }
  }
  OpenSnPETScCall(MatAssemblyBegin(A_shift_, MAT_FINAL_ASSEMBLY));
  OpenSnPETScCall(MatAssemblyEnd(A_shift_, MAT_FINAL_ASSEMBLY));
  OpenSnPETScCall(KSPSetOperators(ksp_, A_shift_, A_shift_));
  const PETScOptionsDatabase::Scope petsc_options_scope(*petsc_options_db_);
  OpenSnPETScCall(KSPSetUp(ksp_));
}

void
CMFDAcceleration::CheckOpticalThickness()
{
  // One-group transport cross section of each coarse cell, weighted by the coarse-group fluxes,
  // times the cell's extent along each axis. The extent is the spread of the face centroids,
  // which is exact for boxes and a fair measure for other cells.
  double local_max[3] = {0.0, 0.0, 0.0};
  std::size_t local_thick = 0;
  for (const auto& coarse_cell : coarse_mesh_.LocalCells())
  {
    const auto& xs = homogenized_xs_.at(coarse_cell.local_id);
    double weighted = 0.0;
    double weight = 0.0;
    for (unsigned int cg = 0; cg < num_coarse_groups_; ++cg)
    {
      const double phi = std::fabs(coarse_phi_new_[coarse_cell.local_id * num_coarse_groups_ + cg]);
      if (xs.diffusion[cg] > 0.0)
      {
        weighted += phi / (3.0 * xs.diffusion[cg]);
        weight += phi;
      }
    }
    if (weight <= 0.0 or coarse_cell.faces.empty())
      continue;
    const double sigma_tr = weighted / weight;
    std::array<double, 3> low = {std::numeric_limits<double>::max(),
                                 std::numeric_limits<double>::max(),
                                 std::numeric_limits<double>::max()};
    std::array<double, 3> high = {std::numeric_limits<double>::lowest(),
                                  std::numeric_limits<double>::lowest(),
                                  std::numeric_limits<double>::lowest()};
    for (const auto& face : coarse_cell.faces)
      for (int d = 0; d < 3; ++d)
      {
        low[d] = std::min(low[d], face.centroid[d]);
        high[d] = std::max(high[d], face.centroid[d]);
      }
    bool thick = false;
    for (int d = 0; d < 3; ++d)
    {
      const double tau = sigma_tr * (high[d] - low[d]);
      local_max[d] = std::max(local_max[d], tau);
      thick = thick or tau > optically_thick_mfp;
    }
    if (thick)
      ++local_thick;
  }
  double global_max[3] = {0.0, 0.0, 0.0};
  opensn::mpi_comm.all_reduce(local_max, 3, global_max, mpi::op::max<double>());
  std::size_t global_thick = 0;
  opensn::mpi_comm.all_reduce(local_thick, global_thick, mpi::op::sum<std::size_t>());
  const int axis = static_cast<int>(std::max_element(global_max, global_max + 3) - global_max);
  const double tau = global_max[axis];
  if (verbose_)
    log.Log() << no_wrap << program_timer.GetTimeString()
              << " CMFD coarse-cell optical thickness: max x = " << global_max[0]
              << ", y = " << global_max[1] << ", z = " << global_max[2] << " mean free paths.";
  if (global_thick == 0)
    return;
  // Convergence rate per update of CMFD with one sweep per update on a 1D slab of coarse cells of
  // thickness tau with scattering ratio 0.99: 0.67 at 1.35, 0.89 at 5.4 and 0.99 at 22 mean free
  // paths. Interpolated in log(tau).
  const double rate =
    tau <= 5.4 ? 0.67 + (0.89 - 0.67) * std::log(tau / 1.35) / std::log(5.4 / 1.35)
               : std::min(0.99, 0.89 + (0.99 - 0.89) * std::log(tau / 5.4) / std::log(22.0 / 5.4));
  const char axis_name = static_cast<char>('x' + axis);
  log.Log0Warning() << no_wrap << "CMFD: " << global_thick << " of "
                    << coarse_mesh_.NumGlobalCells() << " coarse cells are thicker than "
                    << optically_thick_mfp << " transport mean free paths (largest " << tau
                    << " along " << axis_name
                    << "). CMFD with one sweep per update converges slowly on optically thick "
                       "coarse cells, about "
                    << rate
                    << " per update in a scattering medium at this thickness. Refine the "
                       "transport mesh or the coarse mesh along "
                    << axis_name << ".";
}

} // namespace opensn
