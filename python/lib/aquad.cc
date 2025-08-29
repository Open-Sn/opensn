// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/math/quadratures/angular/lebedev_quadrature.h"
#include <pybind11/stl.h>
#include <memory>
#include <stdexcept>

namespace opensn
{

// Wrap quadrature point
void
WrapQuadraturePointPhiTheta(py::module& aquad)
{
  // clang-format off
  py::class_<QuadraturePointPhiTheta> quad_pt_phi_theta(aquad,
    "QuadraturePointPhiTheta",
    R"(
    Angular quadrature point.

    Wrapper of :cpp:class:`opensn::QuadraturePointPhiTheta`.
    )"
  );
  quad_pt_phi_theta.def_readonly(
    "phi",
    &QuadraturePointPhiTheta::phi,
    "Azimuthal angle."
  );
  quad_pt_phi_theta.def_readonly(
    "theta",
    &QuadraturePointPhiTheta::theta,
    "Polar angle."
  );
  quad_pt_phi_theta.def(
    "__repr__",
    [](QuadraturePointPhiTheta& self)
    {
      std::ostringstream os;
      os << "QuadraturePointPhiTheta(phi=" << self.phi << ", theta=" << self.theta << ")";
      return os.str();
    }
  );
  // clang-format on
}

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
  angular_quadrature.def_readonly(
    "abscissae",
    &AngularQuadrature::abscissae,
    "Vector of polar and azimuthal angles."
  );
  angular_quadrature.def_readonly(
    "weights",
    &AngularQuadrature::weights,
    "Quadrature weights."
  );
  angular_quadrature.def_readonly(
    "omegas",
    &AngularQuadrature::omegas,
    "Vector of direction vectors."
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
      [](py::kwargs& params)
      {
        static const std::vector<std::string> required_keys = {"n_polar", "scattering_order"};
        static const std::vector<std::pair<std::string, py::object>> optional_keys = {{"verbose", py::bool_(false)}};
        return construct_from_kwargs<GLProductQuadrature1DSlab, int, int, bool>(params, required_keys, optional_keys);
      }
    ),
    R"(
    Construct a Gauss-Legendre product quadrature for 1D, slab geometry.

    Parameters
    ----------
    n_polar: int
        Number of polar angles.
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
    verbose: bool, default=False
        Verbosity.
    )"
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
      [](py::kwargs& params)
      {
        static const std::vector<std::string> required_keys = {"n_polar", "n_azimuthal", "scattering_order"};
        static const std::vector<std::pair<std::string, py::object>> optional_keys = {{"verbose", py::bool_(false)}};
        return construct_from_kwargs<GLCProductQuadrature2DXY, int, int, int, bool>(params, required_keys, optional_keys);
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
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
    verbose: bool, default=False
        Verbosity.
    )"
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
      [](py::kwargs& params)
      {
        static const std::vector<std::string> required_keys = {"n_polar", "n_azimuthal", "scattering_order"};
        static const std::vector<std::pair<std::string, py::object>> optional_keys = {{"verbose", py::bool_(false)}};
        return construct_from_kwargs<GLCProductQuadrature3DXYZ, int, int, int, bool>(params, required_keys, optional_keys);
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
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
    verbose: bool, default=False
        Verbosity.
    )"
  );
  // clang-format on
}

// Wrap curvilinear product quadrature
void
WrapCurvilinearProductQuadrature(py::module& aquad)
{
  // clang-format off
  // curvilinear product quadrature
  auto curvilinear_product_quadrature = py::class_<CurvilinearProductQuadrature,
                                                   std::shared_ptr<CurvilinearProductQuadrature>,
                                                   ProductQuadrature>(
    aquad,
    "CurvilinearProductQuadrature",
    R"(
    Curvilinear product quadrature.

    Wrapper of :cpp:class:`opensn::CurvilinearProductQuadrature`.
    )"
  );

  // Gauss-Legendre-Chebyshev 2D RZ curvilinear product quadrature
  auto curvilinear_quadrature_glc_2d_rz = py::class_<GLCProductQuadrature2DRZ,
                                                     std::shared_ptr<GLCProductQuadrature2DRZ>,
                                                     CurvilinearProductQuadrature>(
    aquad,
    "GLCProductQuadrature2DRZ",
    R"(
    Gauss-Legendre-Chebyshev product quadrature for 2D, RZ geometry.

    Wrapper of :cpp:class:`opensn::GLCProductQuadrature2DRZ`.
    )"
  );
  curvilinear_quadrature_glc_2d_rz.def(
    py::init(
      [](py::kwargs& params)
      {
        static const std::vector<std::string> required_keys = {"n_polar", "n_azimuthal", "scattering_order"};
        static const std::vector<std::pair<std::string, py::object>> optional_keys = {{"verbose", py::bool_(false)}};
        return construct_from_kwargs<GLCProductQuadrature2DRZ, int, int, int, bool>(params, required_keys, optional_keys);
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
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
    verbose: bool, default=False
        Verbosity.
    )"
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
      [](py::kwargs& params)
      {
        static const std::vector<std::string> required_keys = {"level", "scattering_order"};
        auto [level, scattering_order] = extract_args_tuple<int, int>(params, required_keys);
        std::shared_ptr<SimplifiedLDFESQ::Quadrature> quad(new SimplifiedLDFESQ::Quadrature(scattering_order));
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
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
    )"
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

// Wrap Lebedev quadrature
void
WrapLebedevQuadrature(py::module& aquad)
{
  // clang-format off
  auto lebedev_quadrature = py::class_<LebedevQuadrature,
                                       std::shared_ptr<LebedevQuadrature>,
                                       AngularQuadrature>(
    aquad,
    "LebedevQuadrature",
    R"(
    Lebedev quadrature for spherical integration.
    
    This quadrature provides high-order accuracy for spherical integration with
    symmetric distribution of points on the sphere.

    Wrapper of :cpp:class:`opensn::LebedevQuadrature`.
    )"
  );
  
  lebedev_quadrature.def(
    py::init(
      [](py::kwargs& params)
      {
        static const std::vector<std::string> required_keys = {"order"};
        static const std::vector<std::pair<std::string, py::object>> optional_keys = {{"verbose", py::bool_(false)}};
        return construct_from_kwargs<LebedevQuadrature, int, bool>(params, required_keys, optional_keys);
      }
    ),
    R"(
    Creates a Lebedev quadrature of the specified order.

    Parameters
    ----------
    order: int
        The order of the quadrature.
    verbose: bool, default=False
        Whether to print verbose output during initialization.
    )"
  );
  
  lebedev_quadrature.def(
    "LoadFromOrder",
    &LebedevQuadrature::LoadFromOrder,
    R"(
    Loads quadrature points from an Order.

    Parameters
    ----------
    order: int
        The order of the quadrature.
    verbose: bool, default=False
        Whether to print verbose output during loading.
    )",
    py::arg("order"),
    py::arg("verbose") = false
  );
  // clang-format on
}

// Wrap the angular quadrature components of OpenSn
void
py_aquad(py::module& pyopensn)
{
  py::module aquad = pyopensn.def_submodule("aquad", "Angular quadrature module.");
  WrapQuadraturePointPhiTheta(aquad);
  WrapQuadrature(aquad);
  WrapProductQuadrature(aquad);
  WrapCurvilinearProductQuadrature(aquad);
  WrapSLDFESQuadrature(aquad);
  WrapLebedevQuadrature(aquad);
}

} // namespace opensn
