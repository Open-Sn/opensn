// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"

namespace opensnlua
{

std::shared_ptr<opensn::ProductQuadrature> AQuadCreateGLProductQuadrature1DSlab(int Npolar);

std::shared_ptr<opensn::ProductQuadrature> AQuadCreateGLCProductQuadrature2DXY(int Npolar,
                                                                               int Nazimuthal);

std::shared_ptr<opensn::ProductQuadrature> AQuadCreateGLCProductQuadrature3DXYZ(int Npolar,
                                                                                int Nazimuthal);

std::shared_ptr<opensn::ProductQuadrature> AQuadCreateGLCProductQuadrature2DRZ(int Npolar,
                                                                               int Nazimuthal);

std::shared_ptr<opensn::SimplifiedLDFESQ::Quadrature>
AQuadCreateSLDFESQAngularQuadrature(int level);

void AQuadLocallyRefineSLDFESQ(std::shared_ptr<opensn::SimplifiedLDFESQ::Quadrature> sldfesq,
                               const opensn::Vector3& ref_dir,
                               const double cone_size,
                               const bool dir_as_plane_normal);

void AQuadPrintQuadratureToFile(std::shared_ptr<opensn::SimplifiedLDFESQ::Quadrature> sldfesq,
                                const std::string& file_base);

} // namespace opensnlua
