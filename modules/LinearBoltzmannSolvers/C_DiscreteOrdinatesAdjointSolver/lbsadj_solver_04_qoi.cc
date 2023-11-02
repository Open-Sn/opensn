#include <utility>

#include "lbsadj_solver.h"

#include "mesh/LogicalVolume/LogicalVolume.h"

//###################################################################
/**Subscribes cells to QOIs.*/
size_t lbs::DiscreteOrdinatesAdjointSolver::
  AddResponseFunction(const std::string& qoi_name,
                      std::shared_ptr<chi_mesh::LogicalVolume> logical_volume,
                      const std::string& lua_function_name)
{
  // Make the designation
  ResponseFunctionDesignation qoi_designation(qoi_name,
                                              std::move(logical_volume),
                                              lua_function_name);
  // Make empty subscriber list (will be populated during initialize)
  std::vector<size_t> cell_rf_subscriptions;

  response_functions_.emplace_back(qoi_designation, cell_rf_subscriptions);

  return response_functions_.size() - 1;
}