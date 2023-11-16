#pragma once

namespace opensn
{

template <class PCType, class VecType>
int PreconditionerApplication(PCType pc, VecType vector, VecType action);

} // namespace opensn
