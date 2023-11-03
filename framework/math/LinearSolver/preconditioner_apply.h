#pragma once

namespace chi_math
{
template <class PCType, class VecType>
int PreconditionerApplication(PCType pc, VecType vector, VecType action);
} // namespace chi_math


