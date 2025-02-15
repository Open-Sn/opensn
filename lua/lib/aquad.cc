// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/aquad.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include <memory>

using namespace opensn;

namespace opensnlua
{

std::shared_ptr<ProductQuadrature>
AQuadCreateGLProductQuadrature1DSlab(int Npolar)
{
  bool verbose = false;
  auto quad = std::make_shared<GLProductQuadrature1DSlab>(Npolar, verbose);
  return quad;
}

std::shared_ptr<ProductQuadrature>
AQuadCreateGLCProductQuadrature2DXY(int Npolar, int Nazimuthal)
{
  bool verbose = false;
  auto quad = std::make_shared<GLCProductQuadrature2DXY>(Npolar, Nazimuthal, verbose);
  return quad;
}

std::shared_ptr<ProductQuadrature>
AQuadCreateGLCProductQuadrature3DXYZ(int Npolar, int Nazimuthal)
{
  bool verbose = false;
  auto quad = std::make_shared<GLCProductQuadrature3DXYZ>(Npolar, Nazimuthal, verbose);
  return quad;
}

std::shared_ptr<opensn::ProductQuadrature>
AQuadCreateGLCProductQuadrature2DRZ(int Npolar, int Nazimuthal)
{
  bool verbose = false;
  const auto quad = std::make_shared<GLCProductQuadrature2DRZ>(Npolar, Nazimuthal, verbose);
  return quad;
}

} // namespace opensnlua
