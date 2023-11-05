#pragma once

#include <utility>

#include "opensn/modules/LinearBoltzmannSolvers/B_DiscreteOrdinatesSolver/lbs_discrete_ordinates_solver.h"
#include "opensn/framework/math/chi_math.h"

#include "opensn/modules/LinearBoltzmannSolvers/C_DiscreteOrdinatesAdjointSolver/ResponseFunction/lbs_adj_response_function.h"
#include "opensn/modules/LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

namespace lbs
{

class DiscreteOrdinatesAdjointSolver : public DiscreteOrdinatesSolver
{
protected:
  typedef std::vector<size_t> VecSize_t;
  typedef std::pair<ResponseFunctionDesignation, VecSize_t> RespFuncAndSubs;

public:
  std::vector<std::vector<double>> m_moment_buffers_;

public:
  explicit DiscreteOrdinatesAdjointSolver(const chi::InputParameters& params);
  explicit DiscreteOrdinatesAdjointSolver(const std::string& solver_name);

  DiscreteOrdinatesAdjointSolver(const DiscreteOrdinatesAdjointSolver&) = delete;
  DiscreteOrdinatesAdjointSolver& operator=(const DiscreteOrdinatesAdjointSolver&) = delete;

  void Initialize() override;
  void Execute() override;

  /**
   * Computes the inner product of the flux and the material source.
   */
  double ComputeInnerProduct();

  /**
   * Returns the list of volumetric response functions.
   */
  const std::vector<RespFuncAndSubs>& GetResponseFunctions() const;

  void MakeAdjointXSs();
  void InitQOIs();

  /**
   * Subscribes cells to QOIs.
   */
  size_t AddResponseFunction(const std::string& qoi_name,
                             std::shared_ptr<chi_mesh::LogicalVolume> logical_volume,
                             const std::string& lua_function_name);

  /**
   * Exports an importance map in binary format.
   */
  void ExportImportanceMap(const std::string& file_name);

protected:
  std::vector<RespFuncAndSubs> response_functions_;

public:
  /**
   * Returns the input parameters.
   */
  static chi::InputParameters GetInputParameters();
};

} // namespace lbs
