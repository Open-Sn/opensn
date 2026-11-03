// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/dense_matrix.h"
#include "framework/data_types/vector3.h"
#include <vector>

namespace opensn
{
class Cell;
class CellMapping;

/**
 * Helpers of the blended vacuum-boundary closure ("closure R") of the MIP diffusion operator used
 * by DSA. On a face f with a prescribed incoming flux the closure adds
 * \f[
 *   b_\theta(\delta, v) = \beta_f (1 + \theta_f) \int_f \delta v
 *     - \frac{1 - \theta_f}{2} \int_f D (\partial_n \delta\, v + \partial_n v\, \delta),
 * \f]
 * which is the natural Marshak condition (coefficient \f$2\beta_f \approx 1/2\f$, no gradient
 * terms) for \f$\theta_f = 1\f$ and the half-range transport moment (coefficient \f$\beta_f\f$,
 * half-strength gradient terms) for \f$\theta_f = 0\f$. The blend is
 * \f[
 *   \theta_f = \max\left(\min\left(1, (s_0 / \sigma_t h_f)^2\right), \theta_{coer}\right),
 *   \qquad \frac{(1 - \theta_{coer})^2}{1 + \theta_{coer}} =
 *          \frac{6 \beta_f \sigma_t h_f}{C_{grad}\, n_{f}},
 * \f]
 * where \f$\sigma_t = 1/(3D)\f$, \f$h_f\f$ is the face-dependent orthogonal length, \f$C_{grad}\f$
 * the discrete normal-derivative trace constant of the basis on the face, \f$n_f\f$ the number of
 * faces of the cell carrying the closure and \f$s_0 = 1 + p/2\f$. The second condition is the
 * coercivity bound of the cell form with a safety factor 2, so the operator is SPD for every cell
 * thickness.
 */

/**
 * Half-range moment of a quadrature on a face with unit normal `normal`:
 * \f$\beta_f = \frac12 \sum_d w_d |\Omega_d \cdot n| / \sum_d w_d\f$, the outgoing partial current
 * per unit scalar flux of an isotropic angular flux (1/4 for the exact angular integral).
 */
double QuadratureHalfRangeMoment(const std::vector<Vector3>& omegas,
                                 const std::vector<double>& weights,
                                 const Vector3& normal);

/**
 * Face-dependent orthogonal length \f$h_f\f$ of face `f` of `cell`: \f$d V / A_f\f$ for simplices
 * (d the dimension), \f$V / A_f\f$ for quadrilaterals and 8-vertex hexahedra, \f$2 d V / S\f$ for
 * other polygons and polyhedra (S the total surface), the cell length for slabs.
 */
double FaceOrthogonalLength(const Cell& cell, unsigned int f);

/**
 * Largest generalized eigenvalue \f$\lambda\f$ of \f$N v = \lambda K v\f$ over the range of `K`,
 * with `K` the cell stiffness matrix \f$\int_K \nabla b_i \cdot \nabla b_j\f$ (symmetric positive
 * semi-definite, null space the constants) and `N` the face matrix \f$\int_f \partial_n b_i
 * \partial_n b_j\f$, which also vanishes on constants. Times \f$h_f\f$ this is the trace constant
 * \f$C_{grad}\f$ of the basis on the face.
 */
double MaxEigenvalueOnStiffnessRange(const DenseMatrix<double>& N, const DenseMatrix<double>& K);

/**
 * Trace constant \f$C_{grad} = h_f \max_v \|\partial_n v\|_f^2 / \|\nabla v\|_K^2\f$ of the basis
 * of `cell_mapping` on face `f` of `cell`, from the cell stiffness matrix `stiffness`
 * (\f$\int_K \nabla b_i \cdot \nabla b_j\f$) and the face's normal-derivative mass matrix, with
 * \f$h_f\f$ = FaceOrthogonalLength. Exactly 3 for linear functions on a tetrahedron.
 */
double FaceTraceConstant(const CellMapping& cell_mapping,
                         const Cell& cell,
                         unsigned int f,
                         const DenseMatrix<double>& stiffness);

/**
 * Blend \f$\theta_f\f$ of closure R from \f$\beta_f\f$, \f$\sigma_t h_f\f$, \f$C_{grad}\f$, the
 * number of closure faces of the cell and \f$s_0\f$ (see above). Returns a value in [0, 1].
 */
double ClosureRTheta(
  double beta, double sigma_t_h, double c_grad, unsigned int num_closure_faces, double s0);

} // namespace opensn
