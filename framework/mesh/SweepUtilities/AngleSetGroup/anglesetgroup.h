#pragma once

#include "mesh/SweepUtilities/AngleSet/AngleSet.h"

// ###################################################################
/**Manages the workstages of a single angleset group.*/
class chi_mesh::sweep_management::AngleSetGroup
{
public:
  std::vector<std::shared_ptr<AngleSet>>& AngleSets() { return angle_sets_; }

private:
  std::vector<std::shared_ptr<AngleSet>> angle_sets_;
};
