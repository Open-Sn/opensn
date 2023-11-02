#ifndef CHITECH_LBS_ACCEL_DIFFUSION_H
#define CHITECH_LBS_ACCEL_DIFFUSION_H

#include "acceleration.h"
#include "math/UnknownManager/unknown_manager.h"
#include "petscksp.h"

// ############################################### Forward declarations
namespace chi_mesh
{
class MeshContinuum;
class Cell;
struct Vector3;
} // namespace chi_mesh

namespace chi_math
{
class SpatialDiscretization;
}

namespace lbs
{
struct UnitCellMatrices;
}

namespace lbs::acceleration
{

struct Multigroup_D_and_sigR;

/**Generic diffusion solver for acceleration.*/
class DiffusionSolver
{
protected:
  typedef std::map<int, Multigroup_D_and_sigR> MatID2XSMap;

  const std::string text_name_;
  const chi_mesh::MeshContinuum& grid_;
  const chi_math::SpatialDiscretization& sdm_;
  const chi_math::UnknownManager uk_man_;

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
    double residual_tolerance = 1.0e-4; ///< Residual tol. relative to rhs
    int max_iters = 100;                ///< Maximum iterations
    bool verbose = false;               ///< Verbosity flag
    bool perform_symmetry_check =
      false;                         ///< For debugging only (very expensive)
    std::string source_lua_function; ///< for mms
    std::string ref_solution_lua_function; ///< for mms
    std::string additional_options_string;
    double penalty_factor = 4.0;
  } options;

public:
  DiffusionSolver(std::string text_name,
                  const chi_math::SpatialDiscretization& sdm,
                  const chi_math::UnknownManager& uk_man,
                  std::map<uint64_t, BoundaryCondition> bcs,
                  MatID2XSMap map_mat_id_2_xs,
                  const std::vector<UnitCellMatrices>& unit_cell_matrices,
                  bool verbose,
                  bool requires_ghosts);

  std::string TextName() const;
  const Vec& RHS() const;
  const std::map<uint64_t, BoundaryCondition>& BCS() const {return bcs_;}

  const chi_math::UnknownManager& UnknownStructure() const;
  const chi_math::SpatialDiscretization& SpatialDiscretization() const;

  std::pair<size_t, size_t> GetNumPhiIterativeUnknowns();

  virtual ~DiffusionSolver();

  void Initialize();

  virtual void AssembleAand_b(const std::vector<double>& q_vector) = 0;
  virtual void Assemble_b(const std::vector<double>& q_vector) = 0;
  virtual void Assemble_b(Vec petsc_q_vector) = 0;
  void AddToRHS(const std::vector<double>& values);

  void Solve(std::vector<double>& solution, bool use_initial_guess=false);
  void Solve(Vec petsc_solution, bool use_initial_guess=false);
};

} // namespace lbs::acceleration

#endif // CHITECH_LBS_ACCEL_DIFFUSION_H
