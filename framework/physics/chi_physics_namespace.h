#ifndef CHI_PHYSICS_NAMESPACE_H
#define CHI_PHYSICS_NAMESPACE_H

#include <petscksp.h>

namespace chi_physics
{
enum class OperationType
{
  SINGLE_VALUE = 0,
  FROM_ARRAY = 1,
  SIMPLEXS0 = 20,
  SIMPLEXS1 = 21,
  EXISTING = 22,
  CHI_XSFILE = 23
};

class FieldFunctionGridBased;
class FieldFunctionGridBased;
class Solver;

// 03 Utils
std::string GetPETScConvergedReasonstring(KSPConvergedReason reason);
} // namespace chi_physics

#endif