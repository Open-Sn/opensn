#pragma once

namespace chi_math
{

template <class PCType, class VecType>
struct PreconditionerContext
{
  virtual int PCApply(PCType& pc, VecType& vector, VecType& action) { return 0; }

  virtual ~PreconditionerContext() = default;
};

} // namespace chi_math
