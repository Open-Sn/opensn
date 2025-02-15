// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/product_quadrature.h"

namespace opensnlua
{

std::shared_ptr<opensn::ProductQuadrature> AQuadCreateGLProductQuadrature1DSlab(int Npolar);

std::shared_ptr<opensn::ProductQuadrature> AQuadCreateGLCProductQuadrature2DXY(int Npolar,
                                                                               int Nazimuthal);

std::shared_ptr<opensn::ProductQuadrature> AQuadCreateGLCProductQuadrature3DXYZ(int Npolar,
                                                                                int Nazimuthal);

std::shared_ptr<opensn::ProductQuadrature> AQuadCreateGLCProductQuadrature2DRZ(int Npolar,
                                                                               int Nazimuthal);

} // namespace opensnlua
