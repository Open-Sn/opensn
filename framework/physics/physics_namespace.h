#pragma once

#include <string>
#include <petscksp.h>

namespace opensn
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

/**Gets the string value of a converged reason.*/
std::string GetPETScConvergedReasonstring(KSPConvergedReason reason);

} // namespace opensn
