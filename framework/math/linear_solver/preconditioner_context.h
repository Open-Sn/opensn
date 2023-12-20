#pragma once

namespace opensn
{

template <class PCType, class VecType>
struct PreconditionerContext
{
  virtual int PCApply(PCType& pc, VecType& vector, VecType& action) { return 0; }

  virtual ~PreconditionerContext() = default;
};

} // namespace opensn
