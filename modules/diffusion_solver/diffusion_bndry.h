#pragma once

#include "modules/diffusion_solver/diffusion.h"

namespace chi_diffusion
{
enum class BoundaryType : int
{
  Reflecting = 1,
  Dirichlet = 2,
  Neumann = 3,
  Robin = 4,
  Vacuum = 5
};
}

/**Parent class for diffusion boundaries*/
class chi_diffusion::Boundary
{
public:
  BoundaryType type_;

  explicit Boundary(BoundaryType in_bndry_type) : type_(in_bndry_type) {}
};

/**Reflecting boundary condition.*/
class chi_diffusion::BoundaryReflecting : public chi_diffusion::Boundary
{
public:
  BoundaryReflecting() : Boundary(BoundaryType::Reflecting) {}
};

/**Dirichlet boundary.*/
class chi_diffusion::BoundaryDirichlet : public chi_diffusion::Boundary
{
public:
  double boundary_value = 0.0;

public:
  BoundaryDirichlet() : Boundary(BoundaryType::Dirichlet) {}
  explicit BoundaryDirichlet(double in_bndry_value)
    : Boundary(BoundaryType::Dirichlet), boundary_value(in_bndry_value)
  {
  }
};

/**Robin boundary condition. This type of boundary condition doubles
 * for any boundary condition of the form
 *
\f[
a \phi + b D \hat{n}\cdot \nabla \phi = f
\f]

When \f$ a=0\f$ the boundary condition is equivalent to a <B>Neumann</B>
boundary condition.

\f[
b D\hat{n}\cdot \nabla \phi = f
\f]

 When \f$ a=\frac{1}{4} \f$, \f$ b=\frac{1}{2} \f$
and \f$ f=0 \f$ then the boundary condition is equivalent to a
<B>Vacuum</B> boundary condition.

\f[
\frac{1}{4}\phi + \frac{1}{2}D\hat{n}\cdot \nabla \phi = 0
\f]
 */
class chi_diffusion::BoundaryRobin : public chi_diffusion::Boundary
{
public:
  double a = 0.25;
  double b = 0.5;
  double f = 0.0;

public:
  BoundaryRobin(double a_value, double b_value, double f_value) : Boundary(BoundaryType::Robin)
  {
    a = a_value;
    b = b_value;
    f = f_value;
  }
};
