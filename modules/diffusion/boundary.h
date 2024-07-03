// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <array>

namespace opensn
{
namespace diffusion
{
class Boundary;

enum class BoundaryType : int
{
  Reflecting = 1,
  Dirichlet = 2,
  Neumann = 3,
  Robin = 4,
  Vacuum = 5
};

/**Parent class for diffusion boundaries*/
class Boundary
{
public:
  BoundaryType type_ = BoundaryType::Dirichlet;

  std::array<double, 3> values_ = {0.0, 0.0, 0.0};
};

/** Multigroup diffusion boundary */
class MGBoundary
{
public:
  BoundaryType type_ = BoundaryType::Vacuum;

  std::array<std::vector<double>, 3> mg_values_;
};

} // namespace diffusion
} // namespace opensn

// //###################################################################
// /**Robin boundary condition is a natural (i.w., weak) boundary condition of
// // the form
//  *
// \f[
// a \phi + b D \hat{n}\cdot \nabla \phi = f
// \f]
// When \f$ a=0\f$ the boundary condition is equivalent to a <B>Neumann</B>
// boundary condition.
// \f[
// b D\hat{n}\cdot \nabla \phi = f
// \f]
//  When \f$ a=\frac{1}{4} \f$, \f$ b=\frac{1}{2} \f$
// and \f$ f=0 \f$ then the boundary condition is equivalent to a
// <B>Vacuum</B> boundary condition.
// \f[
// \frac{1}{4}\phi + \frac{1}{2}D\hat{n}\cdot \nabla \phi = 0
// \f]
//  */
// However, one should not set \f$b=0$\f in the hopes of obtaining a
// <B>Dirichlet</B> boundary condition as this type of boundary condition is
// strongly imposed in a different implementation.
