// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"

namespace opensn
{

class TriangleQuadrature : public AngularQuadrature
{
protected:
  unsigned int method_;
  unsigned int sn_;
  unsigned int moments_;

  void TriangleInit();

  void MakeHarmonicIndices(unsigned int scattering_order, int dimension) override;

public:
  explicit TriangleQuadrature(unsigned int method, unsigned int sn, unsigned int moments = 0);

  void BuildDiscreteToMomentOperator(unsigned int scattering_order, int dimension) override;

  void BuildMomentToDiscreteOperator(unsigned int scattering_order, int dimension) override;

  void FilterMoments(unsigned int scattering_order);
};

} // namespace opensn
