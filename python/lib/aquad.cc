// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/logging/log.h"
#include <memory>
#include <stdexcept>

namespace opensn
{

// clang-format off

// Wrap angular quadrature
void wrap_quadrature(py::module &aquad)
{
  // Angular quadrature
  auto angular_quadrature = py::class_<AngularQuadrature, std::shared_ptr<AngularQuadrature>>(
    aquad, "AngularQuadrature",
    R"(
      Angular quadrature.

      Wrapper of :cpp:class:`opensn::AngularQuadrature`.
    )"
  );

  // Product quadrature
  auto product_quadrature = py::class_<ProductQuadrature, std::shared_ptr<ProductQuadrature>, AngularQuadrature>(
    aquad, "ProductQuadrature",
    R"(
      Product quadrature.

      Wrapper of :cpp:class:`opensn::ProductQuadrature`.
    )"
  );

  // Gauss-Legendre 1D slab product quadrature
  auto angular_quadrature_gl_prod_1d_slab = py::class_<GLProductQuadrature1DSlab,
    std::shared_ptr<GLProductQuadrature1DSlab>, ProductQuadrature>(
    aquad, "GLProductQuadrature1DSlab",
    R"(
      Gauss-Legendre product quadrature for 1D, slab geometry.

      Wrapper of :cpp:class:`opensn::GLProductQuadrature1DSlab`.
    )"
  );

  angular_quadrature_gl_prod_1d_slab.def(
    py::init(
      [](int Npolar)
      {
        return std::make_shared<GLProductQuadrature1DSlab>(Npolar);
      }
    ),
    R"(
      Construct a Gauss-Legendre product quadrature for 1D, slab geometry.

      Parameters
      ----------
      Npolar: int
              Number of polar angles.
    )",
    py::arg("Npolar")
  );

  // Gauss-Legendre-Chebyshev 2D XY product quadrature
  auto angular_quadrature_glc_prod_2d_xy = py::class_<GLCProductQuadrature2DXY,
    std::shared_ptr<GLCProductQuadrature2DXY>, ProductQuadrature>(
    aquad, "GLCProductQuadrature2DXY",
    R"(
      Gauss-Legendre-Chebyshev product quadrature for 2D, XY geometry.

      Wrapper of :cpp:class:`opensn::GLCProductQuadrature2DXY`.
    )"
  );

  angular_quadrature_glc_prod_2d_xy.def(
    py::init(
      [](int Npolar, int Nazimuthal)
      {
        return std::make_shared<GLCProductQuadrature2DXY>(Npolar, Nazimuthal);
      }
    ),
    R"(
      Construct a Gauss-Legendre-Chebyshev product quadrature for 2D, XY geometry.

      Parameters
      ----------
      Npolar: int
              Number of polar angles.
      Nazimuthal: int
              Number of azimuthal angles.
    )",
    py::arg("Npolar"), py::arg("Nazimuthal")
  );

  // Gauss-Legendre-Chebyshev 3D XYZ product quadrature
  auto angular_quadrature_glc_prod_3d_xyz = py::class_<GLCProductQuadrature3DXYZ,
    std::shared_ptr<GLCProductQuadrature3DXYZ>, ProductQuadrature>(
    aquad, "GLCProductQuadrature3DXYZ",
    R"(
      Gauss-Legendre-Chebyshev product quadrature for 3D, XYZ geometry.

      Wrapper of :cpp:class:`opensn::GLCProductQuadrature3DXYZ`.
    )"
  );

  angular_quadrature_glc_prod_3d_xyz.def(
    py::init(
      [](int Npolar, int Nazimuthal)
      {
        return std::make_shared<GLCProductQuadrature3DXYZ>(Npolar, Nazimuthal);
      }
    ),
    R"(
      Construct a Gauss-Legendre-Chebyshev product quadrature for 3D, XYZ geometry.

      Parameters
      ----------
      Npolar: int
              Number of polar angles.
      Nazimuthal: int
              Number of azimuthal angles.
    )",
    py::arg("Npolar"), py::arg("Nazimuthal")
  );

  // Simplified LDFESQ quadrature
  auto simplified_LDFESQ_quadrature = py::class_<SimplifiedLDFESQ::Quadrature,
    std::shared_ptr<SimplifiedLDFESQ::Quadrature>, AngularQuadrature>(
    aquad, "SLDFESQuadrature",
    R"(
      Piecewise-linear finite element quadrature using quadrilaterals.

      Wrapper of :cpp:class:`opensn::SimplifiedLDFESQ::Quadrature`.
    )"
  );
}

// Wrap the angular quadrature components of OpenSn
void py_aquad(py::module &pyopensn)
{
  py::module aquad = pyopensn.def_submodule("aquad", "Angular quadrature module.");
  wrap_quadrature(aquad);
}

// clang-format on

} // namespace opensn
