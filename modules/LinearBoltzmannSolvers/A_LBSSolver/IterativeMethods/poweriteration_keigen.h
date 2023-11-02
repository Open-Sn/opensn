#ifndef CHITECH_LBS_POWERITERATION_KEIGEN_H
#define CHITECH_LBS_POWERITERATION_KEIGEN_H

namespace lbs
{
  class LBSSolver;
void PowerIterationKEigen(LBSSolver& lbs_solver,
                          double tolerance,
                          int max_iterations,
                          double& k_eff);
void PowerIterationKEigen1(LBSSolver& lbs_solver,
                           double tolerance,
                           int max_iterations,
                           double& k_eff);
void PowerIterationKEigen2(LBSSolver& lbs_solver,
                           double tolerance,
                           int max_iterations,
                           double& k_eff);

}//namespace lbs

#endif //CHITECH_LBS_POWERITERATION_KEIGEN_H
