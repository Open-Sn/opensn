// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <petscksp.h>

namespace opensn
{

enum class OperationType
{
  SINGLE_VALUE = 0,
  FROM_ARRAY = 1,
  SIMPLE_ONE_GROUP = 20,
  EXISTING = 22,
  OPENSN_XSFILE = 23
};

class FieldFunctionGridBased;
class FieldFunctionGridBased;
class Solver;

/**Gets the string value of a converged reason.*/
std::string GetPETScConvergedReasonstring(KSPConvergedReason reason);

} // namespace opensn
