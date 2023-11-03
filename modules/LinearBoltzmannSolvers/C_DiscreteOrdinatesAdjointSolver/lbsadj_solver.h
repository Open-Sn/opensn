#pragma once

#include <utility>

#include "B_DiscreteOrdinatesSolver/lbs_discrete_ordinates_solver.h"
#include "math/chi_math.h"

#include "ResponseFunction/lbs_adj_response_function.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

namespace lbs
{

// ###################################################################
class DiscreteOrdinatesAdjointSolver : public lbs::DiscreteOrdinatesSolver
{
protected:
  typedef std::vector<size_t> VecSize_t;
  typedef std::pair<ResponseFunctionDesignation, VecSize_t> RespFuncAndSubs;
  std::vector<RespFuncAndSubs> response_functions_;

public:
  std::vector<std::vector<double>> m_moment_buffers_;

public:
  static chi::InputParameters GetInputParameters();
  explicit DiscreteOrdinatesAdjointSolver(const chi::InputParameters& params);
  explicit DiscreteOrdinatesAdjointSolver(const std::string& solver_name);

  DiscreteOrdinatesAdjointSolver(const DiscreteOrdinatesAdjointSolver&) = delete;
  DiscreteOrdinatesAdjointSolver& operator=(const DiscreteOrdinatesAdjointSolver&) = delete;

  double ComputeInnerProduct();
  const std::vector<RespFuncAndSubs>& GetResponseFunctions() const;

  void Initialize() override;
  void MakeAdjointXSs();
  void InitQOIs();
  void Execute() override;

  // 04
  size_t AddResponseFunction(const std::string& qoi_name,
                             std::shared_ptr<chi_mesh::LogicalVolume> logical_volume,
                             const std::string& lua_function_name);
  // 05a
  void ExportImportanceMap(const std::string& file_name);
};

} // namespace lbs

