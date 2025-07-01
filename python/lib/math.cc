// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/math/functions/function.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/vector3.h"
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace opensn
{

// Wrap spherical harmonics
void
WrapYlm(py::module& math)
{
  // clang-format off
  math.def(
    "Ylm",
    &Ylm,
    R"(
    Compute the tesseral spherical harmonics.

    Parameters
    ----------
    l: int
        Degree of the associated Legendre polynomial.
    m: int
        Order of the associated Legendre polynomial.
    theta: float
        Polar angle of the evaluation point.
    varphi: float
        Azimuthal angle of the evaluation point.
    )",
    py::arg("l"),
    py::arg("m"),
    py::arg("theta"),
    py::arg("varphi")
  );
  // clang-format on
}

// Wrap Vector3
void
WrapVector3(py::module& math)
{
  // clang-format off
  auto vector3 = py::class_<Vector3, std::shared_ptr<Vector3>>(
    math,
    "Vector3",
    R"(
    General 3-element vector structure.

    Wrapper of :cpp:class:`opensn::Vector3`.

    Examples
    --------
    >>> a = Vector3(2.0, 3.0, 6.0)
    >>> b = Vector3(2.0, 1.0, -1.0)
    >>> a + b
    Vector3(4, 4, 5)
    >>> a - b
    Vector3(0, 2, 7)
    >>> 2 * a
    Vector3(4, 6, 12)
    >>> a / 2
    Vector3(1, 1.5, 3)
    >>> a *= 2
    >>> a
    Vector3(4, 6, 12)
    >>> a /= 2
    >>> a
    Vector3(2, 3, 6)
    >>> a.Norm()
    7.0
    >>> a @ b  # scalar product
    1.0
    )"
  );
  vector3.def(
    py::init(
      [](double x, double y, double z)
      {
        return std::make_shared<Vector3>(x, y, z);
      }
    ),
    R"(
    Construct a 3-element vector object.

    Parameters
    ----------
    x: float, default=0.0
        X-coordinate.
    y: float, default=0.0
        Y-coordinate.
    z: float, default=0.0
        Z-coordinate.
    )",
    py::arg("x") = 0.0,
    py::arg("y") = 0.0,
    py::arg("z") = 0.0
  );
  vector3.def_readwrite(
    "x",
    &Vector3::x,
    "X-coordinate."
  );
  vector3.def_readwrite(
    "y",
    &Vector3::y,
    "Y-coordinate."
  );
  vector3.def_readwrite(
    "z",
    &Vector3::z,
    "Z-coordinate."
  );
  vector3.def(
    "__add__",
    [](Vector3& self, Vector3& other)
    {
      return self + other;
    }
  );
  vector3.def(
    "__iadd__",
    [](Vector3& self, Vector3& other)
    {
      return self += other;
    }
  );
  vector3.def(
    "__sub__",
    [](Vector3& self, Vector3& other)
    {
      return self - other;
    }
  );
  vector3.def(
    "__isub__",
    [](Vector3& self, Vector3& other)
    {
      return self -= other;
    }
  );
  vector3.def(
    "__mul__",
    [](Vector3& self, double value)
    {
      return self * value;
    }
  );
  vector3.def(
    "__rmul__",
    [](Vector3& self, double value)
    {
      return self * value;
    }
  );
  vector3.def(
    "__imul__",
    [](Vector3& self, double value)
    {
      return self *= value;
    }
  );
  vector3.def(
    "__truediv__",
    [](Vector3& self, double value)
    {
      return self / value;
    }
  );
  vector3.def(
    "__itruediv__",
    [](Vector3& self, double value)
    {
      return self /= value;
    }
  );
  vector3.def(
    "__matmul__",
    [](Vector3& self, Vector3& other)
    {
      return self.Dot(other);
    }
  );
  vector3.def(
    "Norm",
    [](Vector3& self)
    {
      return self.Norm();
    }
  );
  vector3.def(
    "__repr__",
    [](Vector3& self)
    {
      std::ostringstream os;
      os << "Vector3(" << self.x << ", " << self.y << ", " << self.z << ")";
      return os.str();
    }
  );
  // clang-format on
}

// Wrap functors (temporary solution until the Lua interface is eradicated)
void
WrapFunctors(py::module& math)
{
  // clang-format off
  // vector spatial function
  auto vector_spatial_function = py::class_<VectorSpatialFunction, std::shared_ptr<VectorSpatialFunction>>(
    math,
    "VectorSpatialFunction",
    R"(
    Vector spatial function.

    Functions that accept a point and a number of groups as input and return a vector (per group).

    Wrapper of :cpp:class:`opensn::VectorSpatialFunction`.

    Examples
    --------
    >>> # Create from a Python function
    >>> def foo(point, n_groups):
    ...     return [point.x * point.y] * n_groups
    >>> f = VectorSpatialFunction(foo)
    >>> 
    >>> # Create from lambda
    >>> g = VectorSpatialFunction(lambda p, n : [p.x + p.y + p.z] * n)
    >>> 
    >>> # Evaluate
    >>> f(Vector3(1.0, 2.0, 3.0), 2)
    [2.0, 2.0]
    >>> g(Vector3(1.0, 2.0, 3.0), 3)
    [6.0, 6.0, 6.0]
    )"
  );
  vector_spatial_function.def(
    py::init(
      [](const std::function<std::vector<double>(const Vector3&, int)>& func)
      {
        return std::make_shared<VectorSpatialFunction>(func);
      }
    ),
    R"(
    Construct a vector spatial function from associated Python function or lambda.

    Parameters
    ----------
    func: Callable[[pyopensn.math.Vector3, int], List[float]]
        Referenced vector spatial function.
    )",
    py::arg("func")
  );
  vector_spatial_function.def(
    "__call__",
    [](VectorSpatialFunction& self, const Vector3& xyz, int num_groups)
    {
      return self(xyz, num_groups);
    },
    R"(
    Evaluate the associated function.

    Parameters
    ----------
    xyz: pyopensn.math.Vector3
        The xyz coordinates of the point where the function is called.
    num_groups: int
        The number of groups.
    )",
    py::arg("xyz"),
    py::arg("num_groups")
  );
  // clang-format on
}

// Wrap the angular quadrature components of OpenSn
void
py_math(py::module& pyopensn)
{
  py::module math = pyopensn.def_submodule("math", "Math module.");
  WrapYlm(math);
  WrapVector3(math);
  WrapFunctors(math);
}

} // namespace opensn
