#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_adjoint_solver/response_function/lbs_adj_response_function.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/math/functions/spatial_material_function.h"

namespace opensn
{
namespace lbs
{

std::vector<double>
ResponseFunctionDesignation::GetMGResponse(const Cell& cell, const size_t num_groups) const
{

  if (response_function)
    return response_function->Evaluate(cell.centroid_, cell.material_id_, num_groups);
  else
    return std::vector<double>(num_groups, 1.0);
}

} // namespace lbs
} // namespace opensn
