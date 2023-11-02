#ifndef CHITECH_PRECONDITIONER_APPLY_H
#define CHITECH_PRECONDITIONER_APPLY_H

namespace chi_math
{
template <class PCType, class VecType>
int PreconditionerApplication(PCType pc, VecType vector, VecType action);
} // namespace chi_math

#endif // CHITECH_PRECONDITIONER_APPLY_H
