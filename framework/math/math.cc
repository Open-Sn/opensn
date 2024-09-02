// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/math.h"
#include <assert.h>

namespace opensn
{

double
Factorial(const int x)
{
  double factorial_value = 1.0;
  for (int i = 2; i <= x; ++i)
    factorial_value *= i;

  return factorial_value;
}

std::pair<double, double>
OmegaToPhiThetaSafe(const Vector3& omega)
{
  // Notes: asin maps [-1,+1] to [-pi/2,+pi/2]
  //        acos maps [-1,+1] to [0,pi]
  // This mapping requires some logic for determining the azimuthal angle.
  //
  const auto omega_hat = omega.Normalized();

  double mu = omega_hat.z;
  mu = std::min(mu, 1.0);
  mu = std::max(mu, -1.0);

  double theta = acos(mu);

  // Handling omega aligned to k_hat
  if (std::fabs(omega_hat.z) < 1.0e-16)
    return {0.0, theta};

  double cos_phi = omega_hat.x / sin(theta);
  cos_phi = std::min(cos_phi, 1.0);
  cos_phi = std::max(cos_phi, -1.0);

  // Computing varphi for NE and NW quadrant
  if (omega_hat.y >= 0.0)
    return {acos(cos_phi), theta};
  // Computing varphi for SE and SW quadrant
  else
    return {2.0 * M_PI - acos(cos_phi), theta};
}

void
PrintMatrix(const MatDbl& A)
{
  size_t AR = A.size();
  size_t AC = 0;
  if (AR)
    AC = A[0].size();
  else
    std::cout << "A has no rows" << std::endl;

  for (size_t i = 0; i < AR; i++)
  {
    for (size_t j = 0; j < AC; j++)
    {
      std::cout << A[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

void
Scale(MatDbl& A, const double& val)
{
  for (std::vector<double>& Ai : A)
    for (double& Aij : Ai)
      Aij *= val;
}

void
Set(MatDbl& A, const double& val)
{
  for (std::vector<double>& Ai : A)
    for (double& Aij : Ai)
      Aij = val;
}

MatDbl
Transpose(const MatDbl& A)
{
  assert(A.size());
  assert(A[0].size());
  size_t AR = A.size();
  size_t AC = 0;
  if (AR)
    AC = A[0].size();

  MatDbl T(AC, std::vector<double>(AR));
  for (size_t i = 0; i < AR; i++)
    for (size_t j = 0; j < AC; j++)
      T[j][i] = A[i][j];
  return T;
}

void
SwapRow(size_t r1, size_t r2, MatDbl& A)
{
  assert(A.size());
  assert(A[0].size());
  size_t AR = A.size();
  size_t AC = 0;
  if (AR)
    AC = A[0].size();

  assert(r1 >= 0 and r1 < AR and r2 >= 0 and r2 < AR);

  for (size_t j = 0; j < AC; j++)
    std::swap(A[r1][j], A[r2][j]);
}

void
SwapColumn(size_t c1, size_t c2, MatDbl& A)
{
  assert(A.size());
  assert(A[0].size());
  size_t AR = A.size();

  if (A.size())
    assert(c1 >= 0 and c1 < A[0].size() and c2 >= 0 and c2 < A[0].size());

  for (size_t i = 0; i < AR; i++)
    std::swap(A[i][c1], A[i][c2]);
}

MatDbl
MatMul(const MatDbl& A, const double c)
{
  size_t R = A.size();
  size_t C = 0;
  if (R)
    C = A[0].size();

  MatDbl B(R, std::vector<double>(C, 0.));

  for (size_t i = 0; i < R; i++)
    for (size_t j = 0; j < C; j++)
      B[i][j] = A[i][j] * c;

  return B;
}

std::vector<double>
MatMul(const MatDbl& A, const std::vector<double>& x)
{
  size_t R = A.size();
  size_t C = x.size();

  assert(R > 0);
  assert(C == A[0].size());

  std::vector<double> b(R, 0.0);

  for (size_t i = 0; i < R; i++)
  {
    for (size_t j = 0; j < C; j++)
      b[i] += A[i][j] * x[j];
  }

  return b;
}

MatDbl
MatMul(const MatDbl& A, const MatDbl& B)
{
  size_t AR = A.size();

  assert(AR != 0 and B.size() != 0);

  size_t AC = A[0].size();
  size_t BC = B[0].size();

  assert(AC != 0 and BC != 0 and AC == B.size());

  size_t CR = AR;
  size_t CC = BC;
  size_t Cs = AC;

  MatDbl C(CR, std::vector<double>(CC, 0.));

  for (size_t i = 0; i < CR; i++)
    for (size_t j = 0; j < CC; j++)
      for (size_t k = 0; k < Cs; k++)
        C[i][j] += A[i][k] * B[k][j];

  return C;
}

MatDbl
MatAdd(const MatDbl& A, const MatDbl& B)
{
  size_t AR = A.size();
  size_t BR = A.size();

  assert(AR != 0 and B.size() != 0);
  assert(AR == BR);

  size_t AC = A[0].size();
  size_t BC = B[0].size();

  assert(AC != 0 and BC != 0);
  assert(AC == BC);

  MatDbl C(AR, std::vector<double>(AC, 0.0));

  for (size_t i = 0; i < AR; i++)
    for (size_t j = 0; j < AC; j++)
      C[i][j] = A[i][j] + B[i][j];

  return C;
}

MatDbl
MatSubtract(const MatDbl& A, const MatDbl& B)
{
  size_t AR = A.size();
  size_t BR = A.size();

  assert(AR != 0 and B.size() != 0);
  assert(AR == BR);

  size_t AC = A[0].size();
  size_t BC = B[0].size();

  assert(AC != 0 and BC != 0);
  assert(AC == BC);

  MatDbl C(AR, std::vector<double>(AC, 0.0));

  for (size_t i = 0; i < AR; i++)
    for (size_t j = 0; j < AC; j++)
      C[i][j] = B[i][j] - A[i][j];

  return C;
}

double
Determinant(const MatDbl& A)
{
  size_t R = A.size();

  if (R == 1)
    return A[0][0];
  else if (R == 2)
  {
    return A[0][0] * A[1][1] - A[0][1] * A[1][0];
  }
  else if (R == 3)
  {
    return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1] -
           A[0][0] * A[1][2] * A[2][1] - A[0][1] * A[1][0] * A[2][2] - A[0][2] * A[1][1] * A[2][0];
  }
  // http://www.cvl.iis.u-tokyo.ac.jp/~Aiyazaki/tech/teche23.htAl
  else if (R == 4)
  {
    return A[0][0] * A[1][1] * A[2][2] * A[3][3] + A[0][0] * A[1][2] * A[2][3] * A[3][1] +
           A[0][0] * A[1][3] * A[2][1] * A[3][2] + A[0][1] * A[1][0] * A[2][3] * A[3][2] +
           A[0][1] * A[1][2] * A[2][0] * A[3][3] + A[0][1] * A[1][3] * A[2][2] * A[3][0] +
           A[0][2] * A[1][0] * A[2][1] * A[3][3] + A[0][2] * A[1][1] * A[2][3] * A[3][0] +
           A[0][2] * A[1][3] * A[2][0] * A[3][1] + A[0][3] * A[1][0] * A[2][2] * A[3][1] +
           A[0][3] * A[1][1] * A[2][0] * A[3][2] + A[0][3] * A[1][2] * A[2][1] * A[3][0] -
           A[0][0] * A[1][1] * A[2][3] * A[3][2] - A[0][0] * A[1][2] * A[2][1] * A[3][3] -
           A[0][0] * A[1][3] * A[2][2] * A[3][1] - A[0][1] * A[1][0] * A[2][2] * A[3][3] -
           A[0][1] * A[1][2] * A[2][3] * A[3][0] - A[0][1] * A[1][3] * A[2][0] * A[3][2] -
           A[0][2] * A[1][0] * A[2][3] * A[3][1] - A[0][2] * A[1][1] * A[2][0] * A[3][3] -
           A[0][2] * A[1][3] * A[2][1] * A[3][0] - A[0][3] * A[1][0] * A[2][1] * A[3][2] -
           A[0][3] * A[1][1] * A[2][2] * A[3][0] - A[0][3] * A[1][2] * A[2][0] * A[3][1];
  }
  else
  {
    double det = 0;
    for (size_t n = 0; n < R; n++)
    {
      std::vector<std::vector<double>> M = SubMatrix(0, n, A);
      double pm = ((n + 1) % 2) * 2.0 - 1.0;
      det += pm * A[0][n] * Determinant(M);
    }
    return det;
  }
}

MatDbl
SubMatrix(const size_t r, const size_t c, const MatDbl& A)
{
  size_t R = A.size();
  size_t C = 0;
  if (R)
    C = A[0].size();

  assert((r >= 0) and (r < R) and (c >= 0) and (c < C));

  MatDbl B(R - 1, std::vector<double>(C - 1));
  for (size_t i = 0, ii = 0; i < R; ++i)
  {
    if (i != r)
    {
      for (size_t j = 0, jj = 0; j < C; ++j)
      {
        if (j != c)
        {
          B[ii][jj] = A[i][j];
          ++jj;
        }
      }
      ++ii;
    }
  }
  return B;
}

void
GaussElimination(MatDbl& A, std::vector<double>& b, int n)
{
  // Forward elimination
  for (int i = 0; i < n - 1; ++i)
  {
    const std::vector<double>& ai = A[i];
    double bi = b[i];
    double factor = 1.0 / A[i][i];
    for (int j = i + 1; j < n; ++j)
    {
      std::vector<double>& aj = A[j];
      double val = aj[i] * factor;
      b[j] -= val * bi;
      for (int k = i + 1; k < n; ++k)
        aj[k] -= val * ai[k];
    }
  }

  // Back substitution
  for (int i = n - 1; i >= 0; --i)
  {
    const std::vector<double>& ai = A[i];
    double bi = b[i];
    for (int j = i + 1; j < n; ++j)
      bi -= ai[j] * b[j];
    b[i] = bi / ai[i];
  }
}

MatDbl
InverseGEPivoting(const MatDbl& A)
{
  assert(A.size());
  assert(A.size() == A[0].size());

  const unsigned int R = A.size();

  std::vector<std::vector<double>> M(R, std::vector<double>(R, 0.));

  for (unsigned int i = 0; i < R; i++)
    M[i][i] = 1.0;

  std::vector<std::vector<double>> B = A;

  for (unsigned int c = 0; c < R; c++)
  {
    // Find a row with the largest pivot value
    unsigned int max_row = c; // nzr = non-zero row
    for (unsigned int r = c; r < R; ++r)
      if (std::fabs(B[r][c]) > std::fabs(B[max_row][c]))
        max_row = r;

    if (max_row != c)
    {
      SwapRow(max_row, c, B);
      SwapRow(max_row, c, M);
    }

    // Eliminate non-zero values
    for (unsigned int r = 0; r < R; r++)
    {
      if (r != c)
      {
        double g = B[r][c] / B[c][c];
        if (B[r][c] != 0)
        {
          for (unsigned int k = 0; k < R; k++)
          {
            B[r][k] -= B[c][k] * g;
            M[r][k] -= M[c][k] * g;
          }
          B[r][c] = 0;
        }
      }
      else
      {
        double g = 1 / B[c][c];
        for (unsigned int k = 0; k < R; k++)
        {
          B[r][k] *= g;
          M[r][k] *= g;
        }
      }
    }
  }
  return M;
}

MatDbl
Inverse(const MatDbl& A)
{
  size_t R = A.size();
  std::vector<std::vector<double>> M(R, std::vector<double>(R, 0.));
  double f(0.0);

  // Only calculate the determinant if matrix size is less than 5 since
  // the inverse is directly calculated for larger matrices. Otherwise,
  // the inverse routine spends all of its time sitting in the determinant
  // function which is unnecessary.
  if (R < 5)
  {
    f = Determinant(A);
    assert(f != 0.0);
    f = 1.0 / f;
  }

  if (R == 1)
    M[0][0] = f;
  else if (R == 2)
  {
    M[0][0] = A[1][1];
    M[0][1] = -A[0][1];
    M[1][0] = -A[1][0];
    M[1][1] = A[0][0];
    Scale(M, f);
  }
  else if (R == 3)
  {
    M[0][0] = A[2][2] * A[1][1] - A[2][1] * A[1][2];
    M[0][1] = -(A[2][2] * A[0][1] - A[2][1] * A[0][2]);
    M[0][2] = A[1][2] * A[0][1] - A[1][1] * A[0][2];
    M[1][0] = -(A[2][2] * A[1][0] - A[2][0] * A[1][2]);
    M[1][1] = A[2][2] * A[0][0] - A[2][0] * A[0][2];
    M[1][2] = -(A[1][2] * A[0][0] - A[1][0] * A[0][2]);
    M[2][0] = A[2][1] * A[1][0] - A[2][0] * A[1][1];
    M[2][1] = -(A[2][1] * A[0][0] - A[2][0] * A[0][1]);
    M[2][2] = A[1][1] * A[0][0] - A[1][0] * A[0][1];
    Scale(M, f);
  }
  else if (R == 4)
  {
    // http://www.cvl.iis.u-tokyo.ac.jp/~Aiyazaki/tech/teche23.htAl
    M[0][0] = A[1][1] * A[2][2] * A[3][3] + A[1][2] * A[2][3] * A[3][1] +
              A[1][3] * A[2][1] * A[3][2] - A[1][1] * A[2][3] * A[3][2] -
              A[1][2] * A[2][1] * A[3][3] - A[1][3] * A[2][2] * A[3][1];
    M[0][1] = A[0][1] * A[2][3] * A[3][2] + A[0][2] * A[2][1] * A[3][3] +
              A[0][3] * A[2][2] * A[3][1] - A[0][1] * A[2][2] * A[3][3] -
              A[0][2] * A[2][3] * A[3][1] - A[0][3] * A[2][1] * A[3][2];
    M[0][2] = A[0][1] * A[1][2] * A[3][3] + A[0][2] * A[1][3] * A[3][1] +
              A[0][3] * A[1][1] * A[3][2] - A[0][1] * A[1][3] * A[3][2] -
              A[0][2] * A[1][1] * A[3][3] - A[0][3] * A[1][2] * A[3][1];
    M[0][3] = A[0][1] * A[1][3] * A[2][2] + A[0][2] * A[1][1] * A[2][3] +
              A[0][3] * A[1][2] * A[2][1] - A[0][1] * A[1][2] * A[2][3] -
              A[0][2] * A[1][3] * A[2][1] - A[0][3] * A[1][1] * A[2][2];

    M[1][0] = A[1][0] * A[2][3] * A[3][2] + A[1][2] * A[2][0] * A[3][3] +
              A[1][3] * A[2][2] * A[3][0] - A[1][0] * A[2][2] * A[3][3] -
              A[1][2] * A[2][3] * A[3][0] - A[1][3] * A[2][0] * A[3][2];
    M[1][1] = A[0][0] * A[2][2] * A[3][3] + A[0][2] * A[2][3] * A[3][0] +
              A[0][3] * A[2][0] * A[3][2] - A[0][0] * A[2][3] * A[3][2] -
              A[0][2] * A[2][0] * A[3][3] - A[0][3] * A[2][2] * A[3][0];
    M[1][2] = A[0][0] * A[1][3] * A[3][2] + A[0][2] * A[1][0] * A[3][3] +
              A[0][3] * A[1][2] * A[3][0] - A[0][0] * A[1][2] * A[3][3] -
              A[0][2] * A[1][3] * A[3][0] - A[0][3] * A[1][0] * A[3][2];
    M[1][3] = A[0][0] * A[1][2] * A[2][3] + A[0][2] * A[1][3] * A[2][0] +
              A[0][3] * A[1][0] * A[2][2] - A[0][0] * A[1][3] * A[2][2] -
              A[0][2] * A[1][0] * A[2][3] - A[0][3] * A[1][2] * A[2][0];

    M[2][0] = A[1][0] * A[2][1] * A[3][3] + A[1][1] * A[2][3] * A[3][0] +
              A[1][3] * A[2][0] * A[3][1] - A[1][0] * A[2][3] * A[3][1] -
              A[1][1] * A[2][0] * A[3][3] - A[1][3] * A[2][1] * A[3][0];
    M[2][1] = A[0][0] * A[2][3] * A[3][1] + A[0][1] * A[2][0] * A[3][3] +
              A[0][3] * A[2][1] * A[3][0] - A[0][0] * A[2][1] * A[3][3] -
              A[0][1] * A[2][3] * A[3][0] - A[0][3] * A[2][0] * A[3][1];
    M[2][2] = A[0][0] * A[1][1] * A[3][3] + A[0][1] * A[1][3] * A[3][0] +
              A[0][3] * A[1][0] * A[3][1] - A[0][0] * A[1][3] * A[3][1] -
              A[0][1] * A[1][0] * A[3][3] - A[0][3] * A[1][1] * A[3][0];
    M[2][3] = A[0][0] * A[1][3] * A[2][1] + A[0][1] * A[1][0] * A[2][3] +
              A[0][3] * A[1][1] * A[2][0] - A[0][0] * A[1][1] * A[2][3] -
              A[0][1] * A[1][3] * A[2][0] - A[0][3] * A[1][0] * A[2][1];

    M[3][0] = A[1][0] * A[2][2] * A[3][1] + A[1][1] * A[2][0] * A[3][2] +
              A[1][2] * A[2][1] * A[3][0] - A[1][0] * A[2][1] * A[3][2] -
              A[1][1] * A[2][2] * A[3][0] - A[1][2] * A[2][0] * A[3][1];
    M[3][1] = A[0][0] * A[2][1] * A[3][2] + A[0][1] * A[2][2] * A[3][0] +
              A[0][2] * A[2][0] * A[3][1] - A[0][0] * A[2][2] * A[3][1] -
              A[0][1] * A[2][0] * A[3][2] - A[0][2] * A[2][1] * A[3][0];
    M[3][2] = A[0][0] * A[1][2] * A[3][1] + A[0][1] * A[1][0] * A[3][2] +
              A[0][2] * A[1][1] * A[3][0] - A[0][0] * A[1][1] * A[3][2] -
              A[0][1] * A[1][2] * A[3][0] - A[0][2] * A[1][0] * A[3][1];
    M[3][3] = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
              A[0][2] * A[1][0] * A[2][1] - A[0][0] * A[1][2] * A[2][1] -
              A[0][1] * A[1][0] * A[2][2] - A[0][2] * A[1][1] * A[2][0];
    Scale(M, f);
  }
  else
    M = InverseGEPivoting(A);

  return M;
}

double
PowerIteration(const MatDbl& A, std::vector<double>& e_vec, int max_it, double tol)
{
  // Local Variables
  unsigned int n = A.size();
  int it_counter = 0;
  std::vector<double> y(n, 1.0);
  double lambda0 = 0.0;

  // Perform initial iteration outside of loop
  std::vector<double> Ay = MatMul(A, y);
  double lambda = Dot(y, Ay);
  y = VecMul(Ay, 1.0 / Vec2Norm(Ay));
  if (lambda < 0.0)
    Scale(y, -1.0);

  // Perform convergence loop
  bool converged = false;
  while (not converged and it_counter < max_it)
  {
    // Update old eigenvalue
    lambda0 = std::fabs(lambda);
    // Calculate new eigenvalue/eigenvector
    Ay = MatMul(A, y);
    lambda = Dot(y, Ay);
    y = VecMul(Ay, 1.0 / Vec2Norm(Ay));

    // Check if converged or not
    if (std::fabs(std::fabs(lambda) - lambda0) <= tol)
      converged = true;
    // Update counter
    ++it_counter;
  }

  if (lambda < 0.0)
    Scale(y, -1.0);

  // Renormalize eigenvector for the last time
  y = VecMul(Ay, 1.0 / lambda);

  // Set eigenvector, return the eigenvalue
  e_vec = std::move(y);

  return lambda;
}

void
PrintVector(const std::vector<double>& x)
{
  for (auto& xi : x)
    std::cout << xi << ' ';
  std::cout << std::endl;
}

void
Scale(std::vector<double>& x, const double& val)
{
  for (double& xi : x)
    xi *= val;
}

void
Set(std::vector<double>& x, const double& val)
{
  for (double& xi : x)
    xi = val;
}

double
Dot(const std::vector<double>& x, const std::vector<double>& y)
{
  // Error Checks
  assert(x.size() > 0);
  assert(y.size() > 0);
  assert(x.size() == y.size());
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; i++)
    val += x[i] * y[i];

  return val;
}

std::vector<double>
VecMul(const std::vector<double>& x, const double& val)
{
  size_t n = x.size();
  std::vector<double> y(n);

  for (size_t i = 0; i != n; ++i)
    y[i] = val * x[i];

  return y;
}

double
Vec1Norm(const std::vector<double>& x)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; i++)
    val += std::fabs(x[i]);

  return val;
}

double
Vec2Norm(const std::vector<double>& x)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; i++)
    val += x[i] * x[i];

  return std::sqrt(val);
}

double
VecInfinityNorm(const std::vector<double>& x)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; i++)
    val += std::max(std::fabs(x[i]), val);

  return val;
}

double
VecPNorm(const std::vector<double>& x, const double& p)
{
  // Local Variables
  size_t n = x.size();
  double val = 0.0;

  for (size_t i = 0; i != n; i++)
    val += std::pow(std::fabs(x[i]), p);

  return std::pow(val, 1.0 / p);
}

std::vector<double>
operator+(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());
  std::vector<double> result(a.size(), 0.0);

  for (size_t i = 0; i < a.size(); ++i)
    result[i] = a[i] + b[i];

  return result;
}

std::vector<double>
operator-(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());
  std::vector<double> result(a.size(), 0.0);

  for (size_t i = 0; i < a.size(); ++i)
    result[i] = a[i] - b[i];

  return result;
}

double
ComputeL2Change(std::vector<double>& x, std::vector<double>& y)
{
  assert(x.size() == y.size());

  double norm = 0.0;
  for (auto i = 0; i < x.size(); ++i)
  {
    double val = x[i] - y[i];
    norm += val * val;
  }

  double global_norm = 0.0;
  mpi_comm.all_reduce<double>(norm, global_norm, mpi::op::sum<double>());

  return std::sqrt(global_norm);
}

double
ComputePointwiseChange(std::vector<double>& x, std::vector<double>& y)
{
  assert(x.size() == y.size());

  double pw_change = 0.0;
  for (auto i = 0; i < x.size(); ++i)
  {
    double max = std::max(x[i], y[i]);
    double delta = std::fabs(x[i] - y[i]);
    if (max >= std::numeric_limits<double>::min())
      pw_change = std::max(delta / max, pw_change);
    else
      pw_change = std::max(delta, pw_change);
  }

  double global_pw_change = 0.0;
  mpi_comm.all_reduce<double>(pw_change, global_pw_change, mpi::op::max<double>());

  return global_pw_change;
}

} // namespace opensn
