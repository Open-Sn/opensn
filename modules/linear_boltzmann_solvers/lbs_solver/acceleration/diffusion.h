#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/acceleration.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "petscksp.h"

namespace opensn
{
class MeshContinuum;
class Cell;
struct Vector3;
class SpatialDiscretization;

namespace lbs
{
struct UnitCellMatrices;
struct Multigroup_D_and_sigR;

/**
 * Generic diffusion solver for acceleration.
 */
class DiffusionSolver
{
protected:
  typedef std::map<int, Multigroup_D_and_sigR> MatID2XSMap;

  const std::string text_name_;
  const MeshContinuum& grid_;
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

public:
  struct Options
  {
    double residual_tolerance = 1.0e-4;  ///< Residual tol. relative to rhs
    int max_iters = 100;                 ///< Maximum iterations
    bool verbose = false;                ///< Verbosity flag
    bool perform_symmetry_check = false; ///< For debugging only (very expensive)
    std::string additional_options_string;
    double penalty_factor = 4.0;
  } options;

public:
  DiffusionSolver(std::string text_name,
                  const SpatialDiscretization& sdm,
                  const UnknownManager& uk_man,
                  std::map<uint64_t, BoundaryCondition> bcs,
                  MatID2XSMap map_mat_id_2_xs,
                  const std::vector<UnitCellMatrices>& unit_cell_matrices,
                  bool verbose,
                  bool requires_ghosts);

  /**
   * Returns the assigned text name.
   */
  std::string TextName() const;

  /**
   * Returns the right-hand side petsc vector.
   */
  const Vec& RHS() const;

  /**
   * Returns the assigned unknown structure.
   */
  const std::map<uint64_t, BoundaryCondition>& BCS() const { return bcs_; }

  const UnknownManager& UnknownStructure() const;

  /**
   * Returns the associated spatial discretization.
   */
  const class SpatialDiscretization& SpatialDiscretization() const;

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

  /**
   * Adds to the right-hand side without applying spatial discretization.
   */
  void AddToRHS(const std::vector<double>& values);

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

} // namespace lbs
} // namespace opensn
