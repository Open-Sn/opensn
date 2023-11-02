#ifndef CHITECH_SPDS_ADAMSADAMSHAWKINS_H
#define CHITECH_SPDS_ADAMSADAMSHAWKINS_H

#include "SPDS.h"

namespace chi_mesh::sweep_management
{

class SPDS_AdamsAdamsHawkins : public SPDS
{
public:
  SPDS_AdamsAdamsHawkins(const chi_mesh::Vector3& omega,
                         const chi_mesh::MeshContinuum& grid,
                         bool cycle_allowance_flag,
                         bool verbose);
  const std::vector<STDG>& GetGlobalSweepPlanes() const
  {
    return global_sweep_planes_;
  }

private:
  void BuildTaskDependencyGraph(
    const std::vector<std::vector<int>>& global_dependencies,
    bool cycle_allowance_flag);

  std::vector<STDG> global_sweep_planes_; ///< Processor sweep planes
};

}

#endif // CHITECH_SPDS_ADAMSADAMSHAWKINS_H
