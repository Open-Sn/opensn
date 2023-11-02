#include "diffusion_mip.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"

#include <utility>

// ###################################################################
/**Default constructor.*/
lbs::acceleration::DiffusionMIPSolver::DiffusionMIPSolver(
  std::string text_name,
  const chi_math::SpatialDiscretization& sdm,
  const chi_math::UnknownManager& uk_man,
  std::map<uint64_t, BoundaryCondition> bcs,
  MatID2XSMap map_mat_id_2_xs,
  const std::vector<UnitCellMatrices>& unit_cell_matrices,
  const bool verbose /*=false*/)
  : DiffusionSolver(std::move(text_name),
                    sdm,
                    uk_man,
                    std::move(bcs),
                    std::move(map_mat_id_2_xs),
                    unit_cell_matrices,
                    verbose,
                    /*requires_ghosts=*/false)
{
  using SDM_TYPE = chi_math::SpatialDiscretizationType;
  const auto& PWLD = SDM_TYPE ::PIECEWISE_LINEAR_DISCONTINUOUS;

  if (sdm_.Type() != PWLD)
    throw std::logic_error("lbs::acceleration::DiffusionMIPSolver: can only be"
                           " used with PWLD.");
}
