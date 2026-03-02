// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/triangular_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/math/quadratures/angular/lebedev_quadrature.h"
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <limits>
#include <memory>
#include <stdexcept>

namespace opensn
{

// Dictionary for Sn Scattering Source Representation
static std::map<std::string, OperatorConstructionMethod> op_cons_type_map{
  {"standard", OperatorConstructionMethod::STANDARD},
  {"galerkin_one", OperatorConstructionMethod::GALERKIN_ONE},
  {"galerkin_three", OperatorConstructionMethod::GALERKIN_THREE}};

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

// Wrap harmonic indices
static void
WrapHarmonicIndices(py::module& aquad)
{
  py::class_<AngularQuadrature::HarmonicIndices> harmonic_indices(aquad, "HarmonicIndices");
  harmonic_indices.def_readonly("ell", &AngularQuadrature::HarmonicIndices::ell);
  harmonic_indices.def_readonly("m", &AngularQuadrature::HarmonicIndices::m);
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
  angular_quadrature.def(
    "GetDiscreteToMomentOperator",
    [](const AngularQuadrature& self) {
      const auto& op = self.GetDiscreteToMomentOperator();
      if (op.empty()) {
        return py::array_t<double>();
      }
      
      size_t num_rows = op.size();
      size_t num_cols = op[0].size();
      
      // Create numpy array with shape [num_rows, num_cols]
      py::array_t<double> result = py::array_t<double>(
        {num_rows, num_cols},  // shape
        {sizeof(double) * num_cols, sizeof(double)}  // strides (row-major)
      );
      
      py::buffer_info buf = result.request();
      auto* ptr = static_cast<double*>(buf.ptr);
      
      // Copy data row by row
      for (size_t i = 0; i < num_rows; ++i) {
        std::copy(op[i].begin(), op[i].end(), ptr + i * num_cols);
      }
      
      return result;
    },
    "Get the discrete-to-moment operator as a numpy array."
  );
  angular_quadrature.def(
    "GetMomentToDiscreteOperator",
    [](const AngularQuadrature& self) {
      const auto& op = self.GetMomentToDiscreteOperator();
      if (op.empty()) {
        return py::array_t<double>();
      }
      
      size_t num_rows = op.size();
      size_t num_cols = op[0].size();
      
      // Create numpy array with shape [num_rows, num_cols]
      py::array_t<double> result = py::array_t<double>(
        {num_rows, num_cols},  // shape
        {sizeof(double) * num_cols, sizeof(double)}  // strides (row-major)
      );
      
      py::buffer_info buf = result.request();
      auto* ptr = static_cast<double*>(buf.ptr);
      
      // Copy data row by row
      for (size_t i = 0; i < num_rows; ++i) {
        std::copy(op[i].begin(), op[i].end(), ptr + i * num_cols);
      }
      
      return result;
    },
    "Get the moment-to-discrete operator as a numpy array."
  );
  angular_quadrature.def(
    "GetMomentToHarmonicsIndexMap",
    &AngularQuadrature::GetMomentToHarmonicsIndexMap,
    py::return_value_policy::reference_internal
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
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto n_polar = pop_cast(params, "n_polar").cast<unsigned int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        unsigned int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<unsigned int>::max())).cast<unsigned int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<unsigned int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<GLProductQuadrature1DSlab>(n_polar, scattering_order, verbose, method);
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
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
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
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto n_polar = pop_cast(params, "n_polar").cast<unsigned int>();
        auto n_azimuthal = pop_cast(params, "n_azimuthal").cast<unsigned int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        unsigned int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<unsigned int>::max())).cast<unsigned int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<unsigned int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<GLCProductQuadrature2DXY>(n_polar, n_azimuthal, scattering_order, verbose, method);
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
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
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
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto n_polar = pop_cast(params, "n_polar").cast<unsigned int>();
        auto n_azimuthal = pop_cast(params, "n_azimuthal").cast<unsigned int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        unsigned int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<unsigned int>::max())).cast<unsigned int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<unsigned int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<GLCProductQuadrature3DXYZ>(n_polar, n_azimuthal, scattering_order, verbose, method);
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
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
    verbose: bool, default=False
        Verbosity.
    )"
  );
  // clang-format on
}

// Wrap triangular quadrature
void
WrapTriangularQuadrature(py::module& aquad)
{
  // clang-format off
  // triangular quadrature base class
  auto triangular_quadrature = py::class_<TriangularQuadrature, std::shared_ptr<TriangularQuadrature>,
                                          AngularQuadrature>(
    aquad,
    "TriangularQuadrature",
    R"(
    Triangular quadrature base class.

    Unlike product quadratures which have a fixed number of azimuthal angles per polar level,
    triangular quadratures have a varying number of azimuthal angles that decreases
    as the polar angle moves away from the equatorial plane.

    Wrapper of :cpp:class:`opensn::TriangularQuadrature`.
    )"
  );

  // Triangular GLC 3D XYZ quadrature
  auto angular_quadrature_triangular_glc_3d_xyz = py::class_<GLCTriangularQuadrature3DXYZ,
                                                             std::shared_ptr<GLCTriangularQuadrature3DXYZ>,
                                                             TriangularQuadrature>(
    aquad,
    "GLCTriangularQuadrature3DXYZ",
    R"(
    Triangular Gauss-Legendre-Chebyshev quadrature for 3D, XYZ geometry.

    For each polar level away from the equator, there is 1 less azimuthal angle
    per octant. The maximum number of azimuthal angles (at the equator) is
    automatically computed as 2 * n_polar.

    Wrapper of :cpp:class:`opensn::GLCTriangularQuadrature3DXYZ`.
    )"
  );
  angular_quadrature_triangular_glc_3d_xyz.def(
    py::init(
      [](py::kwargs& params)
      {
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto n_polar = pop_cast(params, "n_polar").cast<unsigned int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        unsigned int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<unsigned int>::max())).cast<unsigned int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<unsigned int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<GLCTriangularQuadrature3DXYZ>(n_polar, scattering_order, verbose, method);
      }
    ),
    R"(
    Construct a Triangular Gauss-Legendre-Chebyshev quadrature for 3D, XYZ geometry.

    Parameters
    ----------
    n_polar: int
        Number of polar angles. The maximum number of azimuthal angles (at the equator)
        is automatically computed as ``2 * n_polar``.
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
    verbose: bool, default=False
        Verbosity.
    )"
  );

  // Triangular GLC 2D XY quadrature
  auto angular_quadrature_triangular_glc_2d_xy = py::class_<GLCTriangularQuadrature2DXY,
                                                            std::shared_ptr<GLCTriangularQuadrature2DXY>,
                                                            TriangularQuadrature>(
    aquad,
    "GLCTriangularQuadrature2DXY",
    R"(
    Triangular Gauss-Legendre-Chebyshev quadrature for 2D, XY geometry.

    Only includes points in the upper hemisphere (z >= 0).

    Wrapper of :cpp:class:`opensn::GLCTriangularQuadrature2DXY`.
    )"
  );
  angular_quadrature_triangular_glc_2d_xy.def(
    py::init(
      [](py::kwargs& params)
      {
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto n_polar = pop_cast(params, "n_polar").cast<unsigned int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        unsigned int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<unsigned int>::max())).cast<unsigned int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<unsigned int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<GLCTriangularQuadrature2DXY>(n_polar, scattering_order, verbose, method);
      }
    ),
    R"(
    Construct a Triangular Gauss-Legendre-Chebyshev quadrature for 2D, XY geometry.

    Parameters
    ----------
    n_polar: int
        Number of polar angles (only upper hemisphere will be used). The maximum
        number of azimuthal angles (at the equator) is automatically computed as 2 * n_polar.
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
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
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto n_polar = pop_cast(params, "n_polar").cast<unsigned int>();
        auto n_azimuthal = pop_cast(params, "n_azimuthal").cast<unsigned int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        unsigned int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<unsigned int>::max())).cast<unsigned int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<unsigned int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<GLCProductQuadrature2DRZ>(n_polar, n_azimuthal, scattering_order, verbose, method);
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
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
    verbose: bool, default=False
        Verbosity.
    )"
  );
  // clang-format on
}

// Wrap SLDFES quadrature
void
WrapSLDFEsqQuadrature(py::module& aquad)
{
  // clang-format off
  // Simplified LDFEsq quadrature
  auto sldfesq_quadrature_3d_xyz = py::class_<SLDFEsqQuadrature3DXYZ,
                                              std::shared_ptr<SLDFEsqQuadrature3DXYZ>,
                                              AngularQuadrature>(
    aquad,
    "SLDFEsqQuadrature3DXYZ",
    R"(
    Piecewise-linear finite element quadrature using quadrilaterals.

    Wrapper of :cpp:class:`opensn::SLDFEsqQuadrature3DXYZ`.
    )"
  );
  sldfesq_quadrature_3d_xyz.def(
    py::init(
      [](py::kwargs& params)
      {
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto level = pop_cast(params, "level").cast<int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<int>::max())).cast<int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<SLDFEsqQuadrature3DXYZ>(level, scattering_order, verbose, method);
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
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
    verbose: bool, default=False
        Verbosity.
    )"
  );
  sldfesq_quadrature_3d_xyz.def(
    "LocallyRefine",
    &SLDFEsqQuadrature3DXYZ::LocallyRefine,
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
  sldfesq_quadrature_3d_xyz.def(
    "PrintQuadratureToFile",
    &SLDFEsqQuadrature3DXYZ::PrintQuadratureToFile,
    R"(
    Prints the quadrature to file.

    Parameters
    ----------
    file_base: str
        File base name.
    )",
    py::arg("file_base")
  );

  // 2D SLDFEsq quadrature
  auto sldfesq_quadrature_2d_xy = py::class_<SLDFEsqQuadrature2DXY,
                                             std::shared_ptr<SLDFEsqQuadrature2DXY>,
                                             AngularQuadrature>(
    aquad,
    "SLDFEsqQuadrature2DXY",
    R"(
    Two-dimensional variant of the piecewise-linear finite element quadrature.

    This quadrature is created from the 3D SLDFEsq set by removing directions with negative
    xi.

    Wrapper of :cpp:class:`opensn::SLDFEsqQuadrature2DXY`.
    )"
  );
  sldfesq_quadrature_2d_xy.def(
    py::init(
      [](py::kwargs& params)
      {
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto level = pop_cast(params, "level").cast<int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<int>::max())).cast<int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<SLDFEsqQuadrature2DXY>(level, scattering_order, verbose, method);
      }
    ),
    R"(
    Generates a 2D SLDFEsq quadrature by removing directions with negative xi.

    Parameters
    ----------
    level: int
        Number of subdivisions of the inscribed cube.
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
    verbose: bool, default=False
        Verbosity.
    )"
  );
  sldfesq_quadrature_2d_xy.def(
    "LocallyRefine",
    &SLDFEsqQuadrature2DXY::LocallyRefine,
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
  sldfesq_quadrature_2d_xy.def(
    "PrintQuadratureToFile",
    &SLDFEsqQuadrature2DXY::PrintQuadratureToFile,
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
  // Lebedev 3D XYZ quadrature
  auto angular_quadrature_lebedev_3d_xyz = py::class_<LebedevQuadrature3DXYZ,
                                                     std::shared_ptr<LebedevQuadrature3DXYZ>,
                                                     AngularQuadrature>(
    aquad,
    "LebedevQuadrature3DXYZ",
    R"(
    Lebedev quadrature for 3D, XYZ geometry.

    This quadrature provides high-order accuracy for spherical integration with
    symmetric distribution of points on the sphere.

    Wrapper of :cpp:class:`opensn::LebedevQuadrature3DXYZ`.
    )"
  );

  angular_quadrature_lebedev_3d_xyz.def(
    py::init(
      [](py::kwargs& params)
      {
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto quadrature_order = pop_cast(params, "quadrature_order").cast<unsigned int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        unsigned int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<unsigned int>::max())).cast<unsigned int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<unsigned int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<LebedevQuadrature3DXYZ>(quadrature_order, scattering_order, verbose, method);
      }
    ),
    R"(
    Constructs a Lebedev quadrature for 3D, XYZ geometry.

    Parameters
    ----------
    quadrature_order: int
        The order of the quadrature.
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
    verbose: bool, default=False
        Whether to print verbose output during initialization.
    )"
  );

  // Lebedev 2D XY quadrature
  auto angular_quadrature_lebedev_2d_xy = py::class_<LebedevQuadrature2DXY,
                                                     std::shared_ptr<LebedevQuadrature2DXY>,
                                                     AngularQuadrature>(
    aquad,
    "LebedevQuadrature2DXY",
    R"(
    Lebedev quadrature for 2D, XY geometry.

    This is a 2D version of the Lebedev quadrature that only includes points
    in the upper hemisphere (z >= 0). Points on the equator (z = 0) have their
    weights halved since they are shared between hemispheres.

    Wrapper of :cpp:class:`opensn::LebedevQuadrature2DXY`.
    )"
  );

  angular_quadrature_lebedev_2d_xy.def(
    py::init(
      [](py::kwargs& params)
      {
        auto method_str = pop_cast(params, "operator_method", py::str("standard")).cast<std::string>();
        auto method = op_cons_type_map.at(method_str);
        auto quadrature_order = pop_cast(params, "quadrature_order").cast<unsigned int>();
        auto verbose = pop_cast(params, "verbose", py::bool_(false)).cast<bool>();
        unsigned int scattering_order = 0;
        if (method == OperatorConstructionMethod::GALERKIN_ONE)
          scattering_order = pop_cast(params, "scattering_order",
                                      py::int_(std::numeric_limits<unsigned int>::max())).cast<unsigned int>();
        else
          scattering_order = pop_cast(params, "scattering_order").cast<unsigned int>();
        if (!params.empty())
        {
          std::ostringstream err;
          err << "Unknown argument(s):";
          for (const auto& item : params)
            err << " \"" << py::str(item.first).cast<std::string>() << "\"";
          throw std::runtime_error(err.str());
        }
        return std::make_shared<LebedevQuadrature2DXY>(quadrature_order, scattering_order, verbose, method);
      }
    ),
    R"(
    Constructs a Lebedev quadrature for 2D, XY geometry.

    Parameters
    ----------
    quadrature_order: int
        The order of the quadrature.
    scattering_order: int
        Maximum scattering order supported by the angular quadrature.
        Optional when ``operator_method='galerkin_one'``, in which case the scattering order
        is automatically determined so that the number of moments equals the number of angles.
    operator_method: {'standard', 'galerkin_one', 'galerkin_three'}, default='standard'
        Method used to construct the discrete-to-moment and moment-to-discrete operators.
    verbose: bool, default=False
        Whether to print verbose output during initialization.
    )"
  );
  // clang-format on
}

// Wrap the angular quadrature components of OpenSn
void
py_aquad(py::module& pyopensn)
{
  py::module aquad = pyopensn.def_submodule("aquad", "Angular quadrature module.");
  WrapQuadraturePointPhiTheta(aquad);
  WrapHarmonicIndices(aquad);
  WrapQuadrature(aquad);
  WrapProductQuadrature(aquad);
  WrapTriangularQuadrature(aquad);
  WrapCurvilinearProductQuadrature(aquad);
  WrapSLDFEsqQuadrature(aquad);
  WrapLebedevQuadrature(aquad);
}

} // namespace opensn