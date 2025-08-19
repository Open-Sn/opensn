// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/acceleration.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "petscksp.h"

namespace opensn
{
class MeshContinuum;
class Cell;
struct Vector3;
class SpatialDiscretization;
struct UnitCellMatrices;
struct Multigroup_D_and_sigR;

/// Generic diffusion solver for acceleration.
class DiffusionSolver
{
protected:
  using MatID2XSMap = std::map<int, Multigroup_D_and_sigR>;

  const std::string name_;
  const std::shared_ptr<MeshContinuum> grid_;
  const class SpatialDiscretization& sdm_;
  const UnknownManager uk_man_;

  const std::map<uint64_t, BoundaryCondition> bcs_;

  const MatID2XSMap mat_id_2_xs_map_;

  const std::vector<UnitCellMatrices>& unit_cell_matrices_;

  const int64_t num_local_dofs_;
  const int64_t num_global_dofs_;

  Mat A_ = nullptr;
  Vec rhs_ = nullptr;
  KSP ksp_ = nullptr;

  const bool requires_ghosts_;
  const bool suppress_bcs_;

public:
  struct Options
  {
    /// Residual tol. relative to rhs
    double residual_tolerance = 1.0e-4;
    /// Maximum iterations
    int max_iters = 100;
    /// Verbosity flag
    bool verbose = false;
    /// For debugging only (very expensive)
    bool perform_symmetry_check = false;
    std::string additional_options_string;
    double penalty_factor = 4.0;
  } options;

public:
  DiffusionSolver(std::string name,
                  const SpatialDiscretization& sdm,
                  const UnknownManager& uk_man,
                  std::map<uint64_t, BoundaryCondition> bcs,
                  MatID2XSMap map_mat_id_2_xs,
                  const std::vector<UnitCellMatrices>& unit_cell_matrices,
                  bool suppress_bcs,
                  bool requires_ghosts,
                  bool verbose);

  /// Returns the assigned name.
  std::string GetName() const;

  /// Returns the right-hand side petsc vector.
  const Vec& GetRHS() const;

  /// Returns the assigned unknown structure.
  const std::map<uint64_t, BoundaryCondition>& GetBCS() const { return bcs_; }

  const UnknownManager& GetUnknownStructure() const;

  /// Returns the associated spatial discretization.
  const class SpatialDiscretization& GetSpatialDiscretization() const;

  std::pair<size_t, size_t> GetNumPhiIterativeUnknowns();

  virtual ~DiffusionSolver();

  /**
   * Initializes the diffusion solver. This involves creating the sparse matrix with the appropriate
   * sparsity pattern. Creating the RHS vector. Creating the KSP solver. Setting the very
   * specialized parameters for Hypre's BooomerAMG. Note: `PCSetFromOptions` and `KSPSetFromOptions`
   * are called at the end. Therefore, any number of additional PETSc options can be passed via the
   * commandline.
   */
  void Initialize();

  virtual void AssembleAand_b(const std::vector<double>& q_vector) = 0;
  virtual void Assemble_b(const std::vector<double>& q_vector) = 0;
  virtual void Assemble_b(Vec petsc_q_vector) = 0;

  /// Adds to the right-hand side without applying spatial discretization.
  void AddToRHS(const std::vector<double>& values);

  /// Adds to the entries into the matrix without applying spatial discretization.
  void AddToMatrix(const std::vector<int64_t>& rows,
                   const std::vector<int64_t>& cols,
                   const std::vector<double>& vals);

  /**
   * Solves the system and stores the local solution in the vector provide.
   *
   * \param solution Vector in to which the solution will be parsed.
   * \param use_initial_guess bool [Default:False] Flag, when set, will
   *                 use the values of the output solution as initial guess.
   */
  void Solve(std::vector<double>& solution, bool use_initial_guess = false);

  /**
   * Solves the system and stores the local solution in the vector provide.
   *
   * \param petsc_solution Vector in to which the solution will be parsed.
   * \param use_initial_guess bool [Default:False] Flag, when set, will
   *                 use the values of the output solution as initial guess.
   */
  void Solve(Vec petsc_solution, bool use_initial_guess = false);
};

} // namespace opensn
