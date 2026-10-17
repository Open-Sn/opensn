// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include <string>
#include <vector>

namespace opensn
{
class DiffusionSolver;
class DiscreteOrdinatesProblem;
class LBSGroupset;
class PowerIterationKEigenSolver;
class SpatialDiscretization;
class UnknownManager;
class VectorGhostCommunicator;

/**
 * Base class for LBS acceleration methods.
 */
class DiscreteOrdinatesKEigenAcceleration
{
public:
  explicit DiscreteOrdinatesKEigenAcceleration(const InputParameters& params);

  virtual ~DiscreteOrdinatesKEigenAcceleration() = default;

  /**
   * Public initialize method, to be called by the power iteration solver
   * using this acceleration scheme.
   *
   * Calls the internal Initialize() method.
   */
  void Initialize(PowerIterationKEigenSolver& solver);

  /**
   * Appends one message for each requirement of this acceleration method that the problem
   * violates. The owning solver calls it from its CheckRequirements, so it runs collectively before
   * Initialize() and before each execution. The default has no requirements.
   */
  virtual void CheckRequirements(std::vector<std::string>& errors) const;

  /// The problem this acceleration method was created for.
  const DiscreteOrdinatesProblem& GetProblem() const { return do_problem_; }

  /**
   * Initialization method for derived classes. Must be overridden.
   *
   * Called by the other Initialize() method that takes a solver.
   */
  virtual void Initialize() = 0;

  /**
   * Pre-execute method for derived classes. Must be overridden.
   *
   * Called by the owning solver at the beginning of its Execute().
   */
  virtual void PreExecute() = 0;

  /**
   * Pre-power iteration method for derived classes. Must be overridden.
   *
   * Called by the owning solver at the beginning of the power iteration loop
   * after the fission source is set and the local moments are scaled.
   */
  virtual void PrePowerIteration() = 0;

  /**
   * Post-power iteration method for derived classes. Must be overridden.
   *
   * Called by the owning solver at the end of the power iteration loop
   * after transport is solved.
   *
   * Returns the k-eigenvalue from the acceleration solve.
   */
  virtual double PostPowerIteration() = 0;

  /**
   * Return whether the most recent acceleration update permits the owning power iteration solver
   * to declare convergence. Acceleration schemes can override this when a guarded or rejected
   * update should force at least one additional transport iteration.
   */
  virtual bool AllowsPowerIterationConvergence() const { return true; }

  /**
   * Optional compact convergence information appended to the owning power iteration
   * status line. Acceleration schemes that gate convergence should expose the active
   * gate here so users can tell why PI is still iterating.
   */
  virtual std::string GetPowerIterationConvergenceInfo() const { return {}; }

  const std::string& GetName() const { return name_; }

protected:
  struct GhostInfo
  {
    std::shared_ptr<VectorGhostCommunicator> vector_ghost_communicator = nullptr;
    std::map<uint64_t, uint64_t> ghost_global_id_2_local_map;
  };

  /**
   * Takes an input vector that is the local version of a PWLD discrete space and then
   * makes it continuous by applying nodal averages.
   */
  void NodallyAveragedPWLDVector(const std::vector<double>& input,
                                 std::vector<double>& output) const;

  /// Copies only the scalar moments from an lbs primary flux moments vector.
  void CopyOnlyPhi0(const std::vector<double>& phi_in, std::vector<double>& phi_local);

  /// Copies back only the scalar moments to a lbs primary flux vector.
  void ProjectBackPhi0(const std::vector<double>& input, std::vector<double>& output) const;

  /// Appends an error if an AAH groupset does one sweep per iteration (unstable for this method).
  void RequireMultipleAAHSweeps(std::vector<std::string>& errors) const;

  /// Appends an error if the problem has more than one groupset.
  void RequireSingleGroupset(std::vector<std::string>& errors) const;

  /// The associated DiscreteOrdinatesProblem problem
  DiscreteOrdinatesProblem& do_problem_;

  /// Absolute residual tolerance from parameters
  const double l_abs_tol_;
  /// Maximum allowable iterations from parameters
  const unsigned int max_iters_;
  /// Verbosity flag from parameters
  const bool verbose_;
  /// PETSc options from parameters
  const std::string petsc_options_;
  /// Maximum inner power iteration count, from parameters
  const int pi_max_its_;
  /// k-eigenvalue tolerance for inner power iterations, from parameters
  const double pi_k_tol_;

  /// Groupsets from the LBSProblem
  const std::vector<LBSGroupset>& groupsets_;
  /// Front groupset from the LBSProblem
  LBSGroupset& front_gs_;
  /// Source moments vector from the LBSProblem
  std::vector<double>& q_moments_local_;
  /// Last updated flux vector from the LBSProblem
  std::vector<double>& phi_old_local_;
  /// Newest updated flux vector from the LBSProblem
  std::vector<double>& phi_new_local_;

  /// Associated PI solver, filled during Initialize()
  PowerIterationKEigenSolver* solver_;

  /// Underlying continuous discretization, if needed
  std::shared_ptr<SpatialDiscretization> pwlc_ptr_ = nullptr;
  /// Ghost information, needed with a continuous discretization
  GhostInfo ghost_info_;

  /// Underlying diffusion solver
  std::shared_ptr<DiffusionSolver> diffusion_solver_ = nullptr;

private:
  const std::string name_;

  /**
   * Work vectors
   */
  ///@{
  std::vector<double> copy_only_phi0_tmp_;
  ///@}

public:
  static InputParameters GetInputParameters();

protected:
  static std::vector<uint64_t> MakePWLDGhostIndices(const SpatialDiscretization& pwld,
                                                    const UnknownManager& uk_man);

  static GhostInfo MakePWLDGhostInfo(const SpatialDiscretization& pwld,
                                     const UnknownManager& uk_man);
};

} // namespace opensn
