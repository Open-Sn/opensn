#include "LagrangeDiscontinuous.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Get the number of ghost degrees-of-freedom.*/
size_t LagrangeDiscontinuous::GetNumGhostDOFs(
  const UnknownManager& unknown_manager) const
{
  return 0;
}

// ###################################################################
/**Returns the ghost DOF indices.*/
std::vector<int64_t> LagrangeDiscontinuous::GetGhostDOFIndices(
  const UnknownManager& unknown_manager) const
{
  return {};
}

} // namespace chi_math::spatial_discretization