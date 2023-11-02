#include "PieceWiseLinearContinuous.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Constructor.*/
PieceWiseLinearContinuous::PieceWiseLinearContinuous(const chi_mesh::MeshContinuum& grid,
                                                     chi_math::QuadratureOrder q_order,
                                                     chi_math::CoordinateSystemType cs_type)
  : spatial_discretization::PieceWiseLinearBase(
      grid, q_order, SDMType::PIECEWISE_LINEAR_CONTINUOUS, cs_type)
{
  CreateCellMappings();

  OrderNodes();
}

// ###################################################################
/**Construct a shared object using the protected constructor.*/
std::shared_ptr<PieceWiseLinearContinuous>
PieceWiseLinearContinuous::New(const chi_mesh::MeshContinuum& grid,
                               QuadratureOrder q_order /*=QuadratureOrder::SECOND*/,
                               CoordinateSystemType cs_type /*=CoordinateSystemType::CARTESIAN*/)

{
  const auto PWLC = SpatialDiscretizationType::PIECEWISE_LINEAR_CONTINUOUS;
  // First try to find an existing spatial discretization that matches the
  // one requested.
  for (auto& sdm : Chi::sdm_stack)
    if (sdm->Type() == PWLC and std::addressof(sdm->Grid()) == std::addressof(grid) and
        sdm->GetCoordinateSystemType() == cs_type)
    {
      auto fe_ptr = std::dynamic_pointer_cast<FiniteElementBase>(sdm);

      ChiLogicalErrorIf(not fe_ptr, "Casting failure to FE");

      if (fe_ptr->GetQuadratureOrder() != q_order) break;

      auto sdm_ptr = std::dynamic_pointer_cast<PieceWiseLinearContinuous>(fe_ptr);

      ChiLogicalErrorIf(not sdm_ptr, "Casting failure");

      return sdm_ptr;
    }

  auto new_sdm = std::shared_ptr<PieceWiseLinearContinuous>(
    new PieceWiseLinearContinuous(grid, q_order, cs_type));

  Chi::sdm_stack.push_back(new_sdm);

  return new_sdm;
}

} // namespace chi_math::spatial_discretization
