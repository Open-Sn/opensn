#pragma once

#include "framework/math/quadratures/angular/product_quadrature.h"

namespace opensn
{

/** Base class for curvilinear angular quadratures (product angular
 *  quadratures with additional direction-dependent parameters).
 */
class CurvilinearQuadrature : public ProductQuadrature
{
protected:
  /** Factor to account for angular diamond differencing. */
  std::vector<double> fac_diamond_difference_;

  /** Factor to account for discretisation of the component of the streaming
   *  operator that contains the angular derivative. */
  std::vector<double> fac_streaming_operator_;

  CurvilinearQuadrature() = default;

  virtual ~CurvilinearQuadrature() = default;

public:
  const std::vector<double>& GetDiamondDifferenceFactor() const { return fac_diamond_difference_; }

  const std::vector<double>& GetStreamingOperatorFactor() const { return fac_streaming_operator_; }
};

} // namespace opensn
