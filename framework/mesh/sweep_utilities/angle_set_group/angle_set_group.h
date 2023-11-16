#pragma once

#include "framework/mesh/sweep_utilities/angle_set/angle_set.h"

namespace opensn
{

/**Manages the workstages of a single angleset group.*/
class AngleSetGroup
{
public:
  std::vector<std::shared_ptr<AngleSet>>& AngleSets() { return angle_sets_; }

private:
  std::vector<std::shared_ptr<AngleSet>> angle_sets_;
};

} // namespace opensn
