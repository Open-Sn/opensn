// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/aquad.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/math/quadratures/gausschebyshev_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/math/quadratures/angular/cylindrical_quadrature.h"
#include <cstddef>
#include <memory>

using namespace opensn;

namespace opensnlua
{

std::shared_ptr<ProductQuadrature>
AQuadCreateProductQuadrature(ProductQuadratureType type, int n, int m)
{
  if (type == ProductQuadratureType::GAUSS_LEGENDRE)
  {
    bool verbose = false;
    auto new_quad = std::make_shared<AngularQuadratureProdGL>(n, verbose);
    return new_quad;
  }
  else if (type == ProductQuadratureType::GAUSS_LEGENDRE_CHEBYSHEV)
  {
    bool verbose = false;
    auto new_quad = std::make_shared<AngularQuadratureProdGLC>(n, m, verbose);
    return new_quad;
  }

  opensn::log.LogAllError()
    << "In call to CreateProductQuadrature. Unsupported quadrature type supplied. Given: "
    << (int)type;
  opensn::Exit(EXIT_FAILURE);
  return nullptr;
}

std::shared_ptr<opensn::ProductQuadrature>
AQuadCreateCylindricalProductQuadrature(ProductQuadratureType type, int Np, int Na)
{
  bool verbose = false;
  std::vector<int> vNa;
  vNa.resize(Np, Na);

  switch (type)
  {
    case ProductQuadratureType::GAUSS_LEGENDRE_CHEBYSHEV:
    {
      opensn::log.Log() << "CreateCylindricalProductQuadrature : "
                        << "Creating Gauss-Legendre-Legendre Quadrature\n";

      const auto quad_pol = GaussLegendreQuadrature(Np, verbose);
      std::vector<GaussQuadrature> quad_azi;
      for (const auto& Na : vNa)
        quad_azi.emplace_back(GaussChebyshevQuadrature(Na, verbose));
      const auto new_quad = std::make_shared<CylindricalQuadrature>(quad_pol, quad_azi, verbose);

      return new_quad;
    }

    default:
    {
      opensn::log.LogAllError() << "CreateCylindricalProductQuadrature : "
                                << "Unsupported quadrature type supplied, type="
                                << static_cast<int>(type);
      opensn::Exit(EXIT_FAILURE);
    }
  }
  return nullptr;
}

void
AQuadOptimizeForPolarSymmetry(std::shared_ptr<AngularQuadrature> aquad, double normalization)
{
  if (normalization > 0.0)
    opensn::log.Log() << "Optimizing angular quadrature for polar symmetry. using "
                      << "normalization factor " << normalization << ".";

  aquad->OptimizeForPolarSymmetry(normalization);
}

std::shared_ptr<SimplifiedLDFESQ::Quadrature>
AQuadCreateSLDFESQAngularQuadrature(int level)
{
  auto sldfesq = std::make_shared<SimplifiedLDFESQ::Quadrature>();
  sldfesq->GenerateInitialRefinement(level);
  return sldfesq;
}

void
AQuadLocallyRefineSLDFESQ(std::shared_ptr<SimplifiedLDFESQ::Quadrature> sldfesq,
                          const opensn::Vector3& ref_dir,
                          const double cone_size,
                          const bool dir_as_plane_normal)
{
  sldfesq->LocallyRefine(ref_dir, cone_size, dir_as_plane_normal);
}

void
AQuadPrintQuadratureToFile(std::shared_ptr<SimplifiedLDFESQ::Quadrature> sldfesq,
                           const std::string& file_base)
{
  sldfesq->PrintQuadratureToFile(file_base);
}

} // namespace opensnlua
