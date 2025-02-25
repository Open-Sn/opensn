// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "pyapi.hpp"

#include <memory>
#include <stdexcept>

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/logging/log.h"

namespace opensn {

// Wrap angular qudrature
static void wrap_quadrature(py::module & aquad) {
    // angular quadrature
    auto angular_quadrature = py::class_<AngularQuadrature, std::shared_ptr<AngularQuadrature>>(
        aquad,
        "AngularQuadrature",
        R"(
        Angular quadrature.

        Wrapper of :cpp:class:`opensn::AngularQuadrature`.
        )"
    );
    angular_quadrature.def(
        "OptimizeForPolarSymmetry",
        [](AngularQuadrature & self, double normalization) {
            if (normalization > 0.0) {
                log.Log() << "Optimizing angular quadrature for polar symmetry. using normalization factor "
                          << normalization << ".";
            }
            self.OptimizeForPolarSymmetry(normalization);
        },
        R"(
        Optimizes the angular quadrature for polar symmetry by removing all the direction with downward pointing polar
        angles.

        Parameters
        ----------
        normalization: float, default=-1.0
            A negative number indicates no applied normalization. If a positive number is provided, the weights will be
            normalized to sum to this number.
        )",
        py::arg("normalization") = -1.0
    );
    // product quadrature
    auto product_quadrature = py::class_<ProductQuadrature, std::shared_ptr<ProductQuadrature>, AngularQuadrature>(
        aquad,
        "ProductQuadrature",
        R"(
        Product quadrature.

        Wrapper of :cpp:class:`opensn::ProductQuadrature`.
        )"
    );
    product_quadrature.def(
        py::init(
            [](const std::string & type, int n, int m) {
                std::shared_ptr<ProductQuadrature> result;
                if (type == "gauss-legendre") {
                    result = std::make_shared<AngularQuadratureProdGL>(n, false);
                    return result;
                }
                if (type == "gauss-legendre-chebyshev") {
                    result = std::make_shared<AngularQuadratureProdGLC>(n, m, false);
                    return result;
                }
                log.LogAllError() << "In call to CreateProductQuadrature. Unsupported quadrature type supplied. Given: "
                                  << type;
                throw std::invalid_argument("Unknown option.");
                return std::shared_ptr<ProductQuadrature>(nullptr);
            }
        ),
        R"(
        Construct a product quadrature.

        Parameters
        ----------
        type: {'gauss-legendre', 'gauss-legendre-chebyshev'}
            Quadrature method.
        n: int
            Number of azimuthal angles.
        m: int
            Number of polar angles (for Gauss-Legendre-Chebyshev quadrature).
        )",
        py::arg("type"), py::arg("n"), py::arg("m")
    );
    // Gauss-Legendre product quadrature
    auto angular_quadrature_prod_gl = py::class_<AngularQuadratureProdGL, std::shared_ptr<AngularQuadratureProdGL>, AngularQuadrature>(
        aquad,
        "AngularQuadratureProdGL",
        R"(
        Gauss-Legendre product quadrature.

        Wrapper of :cpp:class:`opensn::AngularQuadratureProdGL`.
        )"
    );
    // Gauss-Legendre-Chebyshev product quadrature
    auto angular_quadrature_prod_glc = py::class_<AngularQuadratureProdGLC, std::shared_ptr<AngularQuadratureProdGLC>, AngularQuadrature>(
        aquad,
        "AngularQuadratureProdGLC",
        R"(
        Gauss-Legendre product quadrature.

        Wrapper of :cpp:class:`opensn::AngularQuadratureProdGLC`.
        )"
    );
    // simplified LDFESQ quadrature
    auto simplified_LDFESQ_quadrature = py::class_<SimplifiedLDFESQ::Quadrature, std::shared_ptr<SimplifiedLDFESQ::Quadrature>, AngularQuadrature>(
        aquad,
        "SLDFESQuadrature",
        R"(
        Piecewise-linear finite element quadrature using quadrilaterals.

        Wrapper of :cpp:class:`opensn::SimplifiedLDFESQ::Quadrature`.
        )"
    );
}

// Wrap the angular quadrature components of OpenSn
void py_aquad(py::module & pyopensn) {
    py::module aquad = pyopensn.def_submodule("aquad", "Angular quadrature module.");
    wrap_quadrature(aquad);
}

}  // namespace opensn
