#pragma once

#include "framework/mesh/sweep_utilities/spds/spds.h"

#include "framework/mesh/sweep_utilities/sweep_namespace.h"

namespace opensn
{
namespace lbs
{

/**Cell-by-Cell (CBC) Sweep Plane Data Structure*/
class CBC_SPDS : public SPDS
{
public:
  CBC_SPDS(const Vector3& omega,
           const MeshContinuum& grid,
           bool cycle_allowance_flag,
           bool verbose);

  const std::vector<Task>& TaskList() const;

protected:
  std::vector<Task> task_list_;
};

} // namespace lbs
} // namespace opensn
