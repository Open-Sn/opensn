#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/angle_set/angle_set.h"

namespace opensn
{
namespace lbs
{

/**Manages the workstages of a single angleset group.*/
class AngleSetGroup
{
private:
  std::vector<std::shared_ptr<AngleSet>> angle_sets_;

public:
  std::vector<std::shared_ptr<AngleSet>>& AngleSets() { return angle_sets_; }
};

} // namespace lbs
} // namespace opensn
