#ifndef CHITECH_LBS_K_RESIDUAL_FUNC_CONTEXT_H
#define CHITECH_LBS_K_RESIDUAL_FUNC_CONTEXT_H

#include <string>
namespace lbs
{
struct KResidualFunctionContext
{
  std::string solver_name;
  double k_eff = 1.0;
};
} // namespace lbs

#endif // CHITECH_LBS_K_RESIDUAL_FUNC_CONTEXT_H
