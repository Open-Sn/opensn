#include "FiniteVolume.h"

#include "chi_runtime.h"
#include "chi_log_exceptions.h"

namespace chi_math::spatial_discretization
{

// ##################################################################
/**Only constructor for this method.*/
FiniteVolume::FiniteVolume(const chi_mesh::MeshContinuum& grid,
                           chi_math::CoordinateSystemType cs_type)
  : SpatialDiscretization(grid, cs_type, SDMType::FINITE_VOLUME)
{
  CreateCellMappings();

  OrderNodes();
}

// ##################################################################
/**Publicly accessible construction handler.*/
std::shared_ptr<FiniteVolume>
FiniteVolume::New(
  const chi_mesh::MeshContinuum& in_grid,
  chi_math::CoordinateSystemType in_cs_type /* = chi_math::CoordinateSystemType::CARTESIAN*/)
{
  // First try to find an existing spatial discretization that matches the
  // one requested.
  for (auto& sdm : Chi::sdm_stack)
    if (sdm->Type() == SpatialDiscretizationType::FINITE_VOLUME and
        std::addressof(sdm->Grid()) == std::addressof(in_grid) and
        sdm->GetCoordinateSystemType() == in_cs_type)
    {
      auto sdm_ptr = std::dynamic_pointer_cast<FiniteVolume>(sdm);

      ChiLogicalErrorIf(not sdm_ptr, "Casting failure");

      return sdm_ptr;
    }

  // If no existing discretization was found then go ahead and make a
  // new one
  auto new_sdm =
    std::shared_ptr<spatial_discretization::FiniteVolume>(new FiniteVolume(in_grid, in_cs_type));

  Chi::sdm_stack.push_back(new_sdm);

  return new_sdm;
}

} // namespace chi_math::spatial_discretization
