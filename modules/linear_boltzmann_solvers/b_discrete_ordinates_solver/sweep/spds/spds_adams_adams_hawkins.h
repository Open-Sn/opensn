#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/spds/spds.h"

namespace opensn
{

class SPDS_AdamsAdamsHawkins : public SPDS
{
public:
  SPDS_AdamsAdamsHawkins(const Vector3& omega,
                         const MeshContinuum& grid,
                         bool cycle_allowance_flag,
                         bool verbose);
  const std::vector<STDG>& GetGlobalSweepPlanes() const { return global_sweep_planes_; }

private:
  /**Builds the task dependency graph.*/
  void BuildTaskDependencyGraph(const std::vector<std::vector<int>>& global_dependencies,
                                bool cycle_allowance_flag);

  std::vector<STDG> global_sweep_planes_; ///< Processor sweep planes
};

} // namespace opensn
