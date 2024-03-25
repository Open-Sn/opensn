#pragma once

#include "framework/math/quadratures/angular/product_quadrature.h"

namespace opensn
{

/** Base class for curvilinear angular quadratures (product angular
 *  quadratures with additional direction-dependent parameters).
 */
class CurvilinearAngularQuadrature : public ProductQuadrature
{
  //  Attributes
protected:
  /** Factor to account for angular diamond differencing. */
  std::vector<double> fac_diamond_difference_;
  /** Factor to account for discretisation of the component of the streaming
   *  operator that contains the angular derivative. */
  std::vector<double> fac_streaming_operator_;

  //  Methods
public:
  /** Return constant reference to fac_diamond_difference. */
  const std::vector<double>& GetDiamondDifferenceFactor() const { return fac_diamond_difference_; }
  /** Return constant reference to fac_streaming_operator. */
  const std::vector<double>& GetStreamingOperatorFactor() const { return fac_streaming_operator_; }

protected:
  /** Default constructor. */
  CurvilinearAngularQuadrature() = default;
  /** Default destructor. */
  virtual ~CurvilinearAngularQuadrature() = default;
};

} // namespace opensn
