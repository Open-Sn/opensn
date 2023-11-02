#include "diffusion_PWLC.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"

namespace lbs::acceleration
{

DiffusionPWLCSolver::DiffusionPWLCSolver(std::string text_name,
                                         const chi_math::SpatialDiscretization& sdm,
                                         const chi_math::UnknownManager& uk_man,
                                         std::map<uint64_t, BoundaryCondition> bcs,
                                         MatID2XSMap map_mat_id_2_xs,
                                         const std::vector<UnitCellMatrices>& unit_cell_matrices,
                                         bool verbose)
  : DiffusionSolver(std::move(text_name),
                    sdm,
                    uk_man,
                    std::move(bcs),
                    std::move(map_mat_id_2_xs),
                    unit_cell_matrices,
                    verbose,
                    /*requires_ghosts=*/true)
{
  using SDM_TYPE = chi_math::SpatialDiscretizationType;
  const auto& PWLC = SDM_TYPE ::PIECEWISE_LINEAR_CONTINUOUS;

  if (sdm_.Type() != PWLC)
    throw std::logic_error("lbs::acceleration::DiffusionPWLCSolver: can only be"
                           " used with PWLC.");
}

} // namespace lbs::acceleration