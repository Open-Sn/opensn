// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include <memory>
#include <stdexcept>

namespace opensn
{

// Wrap angular quadrature
void
WrapQuadrature(py::module& aquad)
{
  // clang-format off
  // angular quadrature
  auto angular_quadrature = py::class_<AngularQuadrature, std::shared_ptr<AngularQuadrature>>(
    aquad,
    "AngularQuadrature",
    R"(
    Angular quadrature.

    Wrapper of :cpp:class:`opensn::AngularQuadrature`.
    )"
  );
  // clang-format on
}

// Wrap product qudrature
void
WrapProductQuadrature(py::module& aquad)
{
  // clang-format off
  // product quadrature
  auto product_quadrature = py::class_<ProductQuadrature, std::shared_ptr<ProductQuadrature>,
                                       AngularQuadrature>(
    aquad,
    "ProductQuadrature",
    R"(
    Product quadrature.

    Wrapper of :cpp:class:`opensn::ProductQuadrature`.
    )"
  );

  // Gauss-Legendre 1D slab product quadrature
  auto angular_quadrature_gl_prod_1d_slab = py::class_<GLProductQuadrature1DSlab,
                                                       std::shared_ptr<GLProductQuadrature1DSlab>,
                                                       ProductQuadrature>(
    aquad,
    "GLProductQuadrature1DSlab",
    R"(
    Gauss-Legendre quadrature for 1D, slab geometry.

    Wrapper of :cpp:class:`opensn::GLProductQuadrature1DSlab`.
    )"
  );
  angular_quadrature_gl_prod_1d_slab.def(
    py::init(
      [](int n_polar, bool verbose)
      {
        return std::make_shared<GLProductQuadrature1DSlab>(n_polar, verbose);
      }
    ),
    R"(
    Construct a Gauss-Legendre product quadrature for 1D, slab geometry.

    Parameters
    ----------
    n_polar: int
        Number of polar angles.
    verbose: bool, default=False
        Verbosity.
    )",
    py::arg("n_polar"),
    py::arg("verbose") = false
  );

  // Gauss-Legendre-Chebyshev 2D XY product quadrature
  auto angular_quadrature_glc_prod_2d_xy = py::class_<GLCProductQuadrature2DXY,
                                                      std::shared_ptr<GLCProductQuadrature2DXY>,
                                                      ProductQuadrature>(
    aquad,
    "GLCProductQuadrature2DXY",
    R"(
    Gauss-Legendre-Chebyshev quadrature for 2D, XY geometry.

    Wrapper of :cpp:class:`opensn::GLCProductQuadrature2DXY`.
    )"
  );
  angular_quadrature_glc_prod_2d_xy.def(
    py::init(
      [](int n_polar, int n_azimuthal, bool verbose)
      {
        return std::make_shared<GLCProductQuadrature2DXY>(n_polar, n_azimuthal, verbose);
      }
    ),
    R"(
    Construct a Gauss-Legendre-Chebyshev product quadrature for 2D, XY geometry.

    Parameters
    ----------
    n_polar: int
        Number of polar angles.
    n_azimuthal: int
        Number of azimuthal angles.
    verbose: bool, default=False
        Verbosity.
    )",
    py::arg("n_polar"),
    py::arg("n_azimuthal"),
    py::arg("verbose") = false
  );

  // Gauss-Legendre-Chebyshev 3D XYZ product quadrature
  auto angular_quadrature_glc_prod_3d_xyz = py::class_<GLCProductQuadrature3DXYZ,
                                                       std::shared_ptr<GLCProductQuadrature3DXYZ>,
                                                       ProductQuadrature>(
    aquad,
    "GLCProductQuadrature3DXYZ",
    R"(
    Gauss-Legendre-Chebyshev quadrature for 3D, XYZ geometry.

    Wrapper of :cpp:class:`opensn::GLCProductQuadrature3DXYZ`.
    )"
  );
  angular_quadrature_glc_prod_3d_xyz.def(
    py::init(
      [](int n_polar, int n_azimuthal)
      {
        return std::make_shared<GLCProductQuadrature3DXYZ>(n_polar, n_azimuthal);
      }
    ),
    R"(
    Construct a Gauss-Legendre-Chebyshev product quadrature for 3D, XYZ geometry.

    Parameters
    ----------
    n_polar: int
        Number of polar angles.
    n_azimuthal: int
        Number of azimuthal angles.
    )",
    py::arg("n_polar"),
    py::arg("n_azimuthal")
  );
  // clang-format on
}

// Wrap curvilinear quadrature
void
WrapCurvilinearQuadrature(py::module& aquad)
{
  // clang-format off
  // curvilinear quadrature
  auto curvilinear_quadrature = py::class_<CurvilinearQuadrature,
                                           std::shared_ptr<CurvilinearQuadrature>,
                                           ProductQuadrature>(
    aquad,
    "CurvilinearQuadrature",
    R"(
    Curvilinear quadrature.

    Wrapper of :cpp:class:`opensn::CurvilinearQuadrature`.
    )"
  );

  // Gauss-Legendre-Chebyshev 2D RZ curvilinear quadrature
  auto curvilinear_quadrature_glc_2d_rz = py::class_<GLCProductQuadrature2DRZ,
                                                     std::shared_ptr<GLCProductQuadrature2DRZ>,
                                                     CurvilinearQuadrature>(
    aquad,
    "GLCProductQuadrature2DRZ",
    R"(
    Gauss-Legendre-Chebyshev quadrature for 2D, RZ geometry.

    Wrapper of :cpp:class:`opensn::GLCProductQuadrature2DRZ`.
    )"
  );
  curvilinear_quadrature_glc_2d_rz.def(
    py::init(
      [](int n_polar, int n_azimuthal, bool verbose)
      {
        return std::make_shared<GLCProductQuadrature2DRZ>(n_polar, n_azimuthal, verbose);
      }
    ),
    R"(
    Construct a Gauss-Legendre Chebyshev product quadrature for 2D, RZ geometry.

    Parameters
    ----------
    n_polar: int
        Number of polar angles.
    n_azimuthal: int
        Number of azimuthal angles.
    verbose: bool, default=False
        Verbosity.
    )",
    py::arg("n_polar"),
    py::arg("n_azimuthal"),
    py::arg("verbose") = false
  );
  // clang-format on
}

// Wrap SLDFES quadrature
void
WrapSLDFESQuadrature(py::module& aquad)
{
  // clang-format off
  // simplified LDFES quadrature
  auto simplified_ldfes_quadrature = py::class_<SimplifiedLDFESQ::Quadrature,
                                                std::shared_ptr<SimplifiedLDFESQ::Quadrature>,
                                                AngularQuadrature>(
    aquad,
    "SLDFESQuadrature",
    R"(
    Piecewise-linear finite element quadrature using quadrilaterals.

    Wrapper of :cpp:class:`opensn::SimplifiedLDFESQ::Quadrature`.
    )"
  );
  simplified_ldfes_quadrature.def(
    py::init(
      [](int level)
      {
        SimplifiedLDFESQ::Quadrature * quad = new SimplifiedLDFESQ::Quadrature();
        quad->GenerateInitialRefinement(level);
        return quad;
      }
    ),
    R"(
    Generates uniform spherical quadrilaterals from the subdivision of an inscribed cube.

    Parameters
    ----------
    level: int
        Number of subdivisions of the inscribed cube.
    )",
    py::arg("level")
  );
  simplified_ldfes_quadrature.def(
    "LocallyRefine",
    &SimplifiedLDFESQ::Quadrature::LocallyRefine,
    R"(
    Locally refines the cells.

    Parameters
    ----------
    ref_dir: pyopensn.math.Vector3
        Reference direction :math:`\vec{r}`.
    cone_size: float
        Cone size (in radians) :math:`\theta`.
    dir_as_plane_normal: bool, default=False
        If true, interpret SQ-splitting as when :math:`|\omega \cdot \vec{r}| < \sin(\theta)`.
        Otherwise, SQs will be split if :math:`\omega \cdot \vec{r} > \cos(\theta)`.
    )",
    py::arg("ref_dir"),
    py::arg("cone_size"),
    py::arg("dir_as_plane_normal") = false
  );
  simplified_ldfes_quadrature.def(
    "PrintQuadratureToFile",
    &SimplifiedLDFESQ::Quadrature::PrintQuadratureToFile,
    R"(
    Prints the quadrature to file.

    Parameters
    ----------
    file_base: str
        File base name.
    )",
    py::arg("file_base")
  );
  // clang-format on
}

// Wrap the angular quadrature components of OpenSn
void
py_aquad(py::module& pyopensn)
{
  py::module aquad = pyopensn.def_submodule("aquad", "Angular quadrature module.");
  WrapQuadrature(aquad);
  WrapProductQuadrature(aquad);
  WrapCurvilinearQuadrature(aquad);
  WrapSLDFESQuadrature(aquad);
}

} // namespace opensn
