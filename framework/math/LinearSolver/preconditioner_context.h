#ifndef CHITECH_PRECONDITIONER_CONTEXT_H
#define CHITECH_PRECONDITIONER_CONTEXT_H

namespace chi_math
{

template<class PCType, class VecType>
struct PreconditionerContext
{
  virtual int PCApply(PCType& pc, VecType& vector, VecType& action)
  {return 0;}

  virtual ~PreconditionerContext() = default;
};

}//namespace chi_math

#endif //CHITECH_PRECONDITIONER_CONTEXT_H
