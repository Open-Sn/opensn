// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/aquad.h"
#include <memory>

using namespace opensn;

namespace opensnlua
{

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
