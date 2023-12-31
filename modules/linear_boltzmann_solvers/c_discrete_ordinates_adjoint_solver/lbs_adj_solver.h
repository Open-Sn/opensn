#pragma once

#include <utility>

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "framework/math/math.h"

#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_adjoint_solver/response_function/lbs_adj_response_function.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"

namespace opensn
{
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
  explicit DiscreteOrdinatesAdjointSolver(const InputParameters& params);
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
                             std::shared_ptr<LogicalVolume> logical_volume,
                             std::shared_ptr<ResponseFunction> function);

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
  static InputParameters GetInputParameters();
};

} // namespace lbs
} // namespace opensn
