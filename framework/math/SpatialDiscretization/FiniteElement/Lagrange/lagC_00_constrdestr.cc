#include "LagrangeContinuous.h"

#include "chi_log_exceptions.h"

namespace chi_math::spatial_discretization
{

LagrangeContinuous::LagrangeContinuous(const chi_mesh::MeshContinuum& grid,
                                       QuadratureOrder q_order,
                                       CoordinateSystemType cs_type)
  : LagrangeBase(grid, q_order, SpatialDiscretizationType::LAGRANGE_CONTINUOUS, cs_type)
{
  OrderNodes();
}

// ###################################################################
/**Construct a shared object using the protected constructor.*/
std::shared_ptr<LagrangeContinuous>
LagrangeContinuous::New(const chi_mesh::MeshContinuum& grid,
                        QuadratureOrder q_order /*=QuadratureOrder::SECOND*/,
                        CoordinateSystemType cs_type /*=CoordinateSystemType::CARTESIAN*/)

{
  const auto LagrangeC = SpatialDiscretizationType::LAGRANGE_CONTINUOUS;
  // First try to find an existing spatial discretization that matches the
  // one requested.
  for (auto& sdm : Chi::sdm_stack)
    if (sdm->Type() == LagrangeC and std::addressof(sdm->Grid()) == std::addressof(grid) and
        sdm->GetCoordinateSystemType() == cs_type)
    {
      auto fe_ptr = std::dynamic_pointer_cast<FiniteElementBase>(sdm);

      ChiLogicalErrorIf(not fe_ptr, "Casting failure to FE");

      if (fe_ptr->GetQuadratureOrder() != q_order) break;

      auto sdm_ptr = std::dynamic_pointer_cast<LagrangeContinuous>(fe_ptr);

      ChiLogicalErrorIf(not sdm_ptr, "Casting failure");

      return sdm_ptr;
    }

  auto new_sdm =
    std::shared_ptr<LagrangeContinuous>(new LagrangeContinuous(grid, q_order, cs_type));

  Chi::sdm_stack.push_back(new_sdm);

  return new_sdm;
}

} // namespace chi_math::spatial_discretization