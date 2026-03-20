Introduction to Angular Quadratures
===================================

Discrete-ordinates codes need angular quadratures to select the sweeping directions and to compute moments of the angular flux.

Several choices are available:

- A Product Quadrature set (using Gauss-Legendre quadrature along the polar angle, and a Gauss-Chebyshev quadrature along the azimuthal angle)
- A Triangular Quadrature set (using Gauss-Legendre quadrature along the polar angle, with decreasing orders of Gauss-Chebyshev quadrature along the azimuthal angle)
- A Linear Discontinuous Finite Element (LDFE) Quadrature set that allows for local angular refinement
- A Lebedev Quadrature set


Additionally, there are various options for how the discrete-to-moment flux map is constructed. 

- The Standard construction method
- Galerkin-Quadrature Method 1
- Galerkin-Quadrature Method 3
