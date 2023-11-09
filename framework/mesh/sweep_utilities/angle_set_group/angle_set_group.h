#pragma once

#include "framework/mesh/sweep_utilities/angle_set/angle_set.h"

/**Manages the workstages of a single angleset group.*/
class chi_mesh::sweep_management::AngleSetGroup
{
public:
  std::vector<std::shared_ptr<AngleSet>>& AngleSets() { return angle_sets_; }

private:
  std::vector<std::shared_ptr<AngleSet>> angle_sets_;
};
