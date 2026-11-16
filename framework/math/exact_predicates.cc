// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/exact_predicates.h"
#include <cmath>
#include <limits>
#include <array>

namespace opensn
{

namespace
{

/**
 * A floating-point expansion: a sum of doubles that are nonoverlapping and sorted by increasing
 * magnitude, so that the sign of the expansion is the sign of its last nonzero component. The
 * operations are the error-free transformations of Shewchuk (1997); products use a fused
 * multiply-add for the exact rounding error.
 */
/// At most 192 components are needed for the exact orientation determinant (Shewchuk 1997).
struct Expansion
{
  // Left uninitialized on purpose: only the first n components are ever read, and the expansion
  // is created on every exact evaluation.
  std::array<double, 256> c; // NOLINT(cppcoreguidelines-pro-type-member-init)
  int n = 0;
  double operator[](const int k) const { return c[k]; }
  bool empty() const { return n == 0; }
  int size() const { return n; }
  void push_back(const double v) { c[n++] = v; }
};

inline void
TwoSum(const double a, const double b, double& x, double& y)
{
  x = a + b;
  const double bvirt = x - a;
  const double avirt = x - bvirt;
  const double bround = b - bvirt;
  const double around = a - avirt;
  y = around + bround;
}

/// TwoSum for |a| >= |b|.
inline void
FastTwoSum(const double a, const double b, double& x, double& y)
{
  x = a + b;
  const double bvirt = x - a;
  y = b - bvirt;
}

inline void
TwoDiff(const double a, const double b, double& x, double& y)
{
  x = a - b;
  const double bvirt = a - x;
  const double avirt = x + bvirt;
  const double bround = bvirt - b;
  const double around = a - avirt;
  y = around + bround;
}

inline void
TwoProduct(const double a, const double b, double& x, double& y)
{
  x = a * b;
  y = std::fma(a, b, -x);
}

void
Difference(const double a, const double b, Expansion& e)
{
  double x = 0.0;
  double y = 0.0;
  TwoDiff(a, b, x, y);
  e.n = 0;
  e.push_back(y);
  e.push_back(x);
}

/// Sum of two expansions (Shewchuk's fast expansion sum with zero elimination).
void
Sum(const Expansion& e, const Expansion& f, Expansion& h)
{
  h.n = 0;
  if (e.empty())
  {
    h = f;
    return;
  }
  if (f.empty())
  {
    h = e;
    return;
  }
  int eindex = 0;
  int findex = 0;
  double enow = e[0];
  double fnow = f[0];
  double Q = 0.0;
  // Take the component of smaller magnitude first.
  const auto take_e = [&]() { return (fnow > enow) == (fnow > -enow); };
  if (take_e())
  {
    Q = enow;
    ++eindex;
  }
  else
  {
    Q = fnow;
    ++findex;
  }
  double Qnew = 0.0;
  double hh = 0.0;
  if (eindex < e.size() and findex < f.size())
  {
    enow = e[eindex];
    fnow = f[findex];
    if (take_e())
    {
      FastTwoSum(enow, Q, Qnew, hh);
      ++eindex;
    }
    else
    {
      FastTwoSum(fnow, Q, Qnew, hh);
      ++findex;
    }
    Q = Qnew;
    if (hh != 0.0)
      h.push_back(hh);
    while (eindex < e.size() and findex < f.size())
    {
      enow = e[eindex];
      fnow = f[findex];
      if (take_e())
      {
        TwoSum(Q, enow, Qnew, hh);
        ++eindex;
      }
      else
      {
        TwoSum(Q, fnow, Qnew, hh);
        ++findex;
      }
      Q = Qnew;
      if (hh != 0.0)
        h.push_back(hh);
    }
  }
  while (eindex < e.size())
  {
    TwoSum(Q, e[eindex], Qnew, hh);
    ++eindex;
    Q = Qnew;
    if (hh != 0.0)
      h.push_back(hh);
  }
  while (findex < f.size())
  {
    TwoSum(Q, f[findex], Qnew, hh);
    ++findex;
    Q = Qnew;
    if (hh != 0.0)
      h.push_back(hh);
  }
  if (Q != 0.0 or h.empty())
    h.push_back(Q);
}

/// Product of an expansion and a double (Shewchuk's scale expansion with zero elimination).
void
Scale(const Expansion& e, const double b, Expansion& h)
{
  h.n = 0;
  if (e.empty())
    return;
  double Q = 0.0;
  double hh = 0.0;
  TwoProduct(e[0], b, Q, hh);
  if (hh != 0.0)
    h.push_back(hh);
  for (int i = 1; i < e.size(); ++i)
  {
    double product1 = 0.0;
    double product0 = 0.0;
    TwoProduct(e[i], b, product1, product0);
    double sum = 0.0;
    TwoSum(Q, product0, sum, hh);
    if (hh != 0.0)
      h.push_back(hh);
    FastTwoSum(product1, sum, Q, hh);
    if (hh != 0.0)
      h.push_back(hh);
  }
  if (Q != 0.0 or h.empty())
    h.push_back(Q);
}

/// Product of two expansions.
void
Product(const Expansion& e, const Expansion& f, Expansion& result)
{
  result.n = 0;
  Expansion scaled;
  Expansion partial;
  for (int k = 0; k < f.size(); ++k)
  {
    Scale(e, f[k], scaled);
    Sum(result, scaled, partial);
    result = partial;
  }
}

void
Negate(Expansion& e)
{
  for (int k = 0; k < e.size(); ++k)
    e.c[k] = -e.c[k];
}

int
Sign(const Expansion& e)
{
  for (int k = e.size() - 1; k >= 0; --k)
  {
    if (e[k] > 0.0)
      return 1;
    if (e[k] < 0.0)
      return -1;
  }
  return 0;
}

int
Orient3DExact(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d)
{
  Expansion adx, ady, adz, bdx, bdy, bdz, cdx, cdy, cdz;
  Difference(a.x, d.x, adx);
  Difference(a.y, d.y, ady);
  Difference(a.z, d.z, adz);
  Difference(b.x, d.x, bdx);
  Difference(b.y, d.y, bdy);
  Difference(b.z, d.z, bdz);
  Difference(c.x, d.x, cdx);
  Difference(c.y, d.y, cdy);
  Difference(c.z, d.z, cdz);
  Expansion u, v, m1, m2, m3, t1, t2, t3, s12, det;
  Product(bdx, cdy, u);
  Product(cdx, bdy, v);
  Negate(v);
  Sum(u, v, m1);
  Product(cdx, ady, u);
  Product(adx, cdy, v);
  Negate(v);
  Sum(u, v, m2);
  Product(adx, bdy, u);
  Product(bdx, ady, v);
  Negate(v);
  Sum(u, v, m3);
  Product(adz, m1, t1);
  Product(bdz, m2, t2);
  Product(cdz, m3, t3);
  Sum(t1, t2, s12);
  Sum(s12, t3, det);
  return Sign(det);
}

inline void
TwoDiffTail(const double a, const double b, const double x, double& y)
{
  const double bvirt = a - x;
  const double avirt = x + bvirt;
  const double bround = bvirt - b;
  const double around = a - avirt;
  y = around + bround;
}

/// The exact difference of two exact products of doubles as a four-component expansion.
void
TwoTwoProductDiff(const double a, const double b, const double c, const double d, Expansion& out)
{
  double ab1 = 0.0, ab0 = 0.0, cd1 = 0.0, cd0 = 0.0;
  TwoProduct(a, b, ab1, ab0);
  TwoProduct(c, d, cd1, cd0);
  // Two_Two_Diff(ab1, ab0, cd1, cd0): x0 .. x3 in increasing magnitude.
  double i = 0.0, x0 = 0.0, j = 0.0, zero = 0.0, x1 = 0.0, x3 = 0.0, x2 = 0.0;
  TwoDiff(ab0, cd0, i, x0);
  TwoSum(ab1, i, j, zero);
  TwoDiff(zero, cd1, i, x1);
  TwoSum(j, i, x3, x2);
  out.n = 0;
  out.push_back(x0);
  out.push_back(x1);
  out.push_back(x2);
  out.push_back(x3);
}

/**
 * Stages B and C of Shewchuk's adaptive orient3d, then the exact evaluation: the determinant of
 * the rounded differences is computed exactly and corrected with the round-off of the
 * differences, each stage with a rigorous error bound.
 */
int
Orient3DAdaptive(
  const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d, const double permanent)
{
  constexpr double epsilon = std::numeric_limits<double>::epsilon() / 2.0;
  constexpr double resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
  constexpr double errboundB = (3.0 + 28.0 * epsilon) * epsilon;
  constexpr double errboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;

  const double adx = a.x - d.x;
  const double bdx = b.x - d.x;
  const double cdx = c.x - d.x;
  const double ady = a.y - d.y;
  const double bdy = b.y - d.y;
  const double cdy = c.y - d.y;
  const double adz = a.z - d.z;
  const double bdz = b.z - d.z;
  const double cdz = c.z - d.z;

  // Stage B: the determinant of the rounded differences, exactly.
  Expansion bc, ca, ab, adet, bdet, cdet, abdet, fin;
  TwoTwoProductDiff(bdx, cdy, cdx, bdy, bc);
  TwoTwoProductDiff(cdx, ady, adx, cdy, ca);
  TwoTwoProductDiff(adx, bdy, bdx, ady, ab);
  Scale(bc, adz, adet);
  Scale(ca, bdz, bdet);
  Scale(ab, cdz, cdet);
  Sum(adet, bdet, abdet);
  Sum(abdet, cdet, fin);
  double det = 0.0;
  for (int k = 0; k < fin.size(); ++k)
    det += fin[k];
  double errbound = errboundB * permanent;
  if (det >= errbound or -det >= errbound)
    return det > 0.0 ? 1 : det < 0.0 ? -1 : 0;

  // Stage C: first-order correction with the round-off of the differences.
  double adxtail = 0.0, bdxtail = 0.0, cdxtail = 0.0, adytail = 0.0, bdytail = 0.0, cdytail = 0.0,
         adztail = 0.0, bdztail = 0.0, cdztail = 0.0;
  TwoDiffTail(a.x, d.x, adx, adxtail);
  TwoDiffTail(b.x, d.x, bdx, bdxtail);
  TwoDiffTail(c.x, d.x, cdx, cdxtail);
  TwoDiffTail(a.y, d.y, ady, adytail);
  TwoDiffTail(b.y, d.y, bdy, bdytail);
  TwoDiffTail(c.y, d.y, cdy, cdytail);
  TwoDiffTail(a.z, d.z, adz, adztail);
  TwoDiffTail(b.z, d.z, bdz, bdztail);
  TwoDiffTail(c.z, d.z, cdz, cdztail);
  if (adxtail == 0.0 and bdxtail == 0.0 and cdxtail == 0.0 and adytail == 0.0 and bdytail == 0.0 and
      cdytail == 0.0 and adztail == 0.0 and bdztail == 0.0 and cdztail == 0.0)
    return det > 0.0 ? 1 : det < 0.0 ? -1 : 0;
  errbound = errboundC * permanent + resulterrbound * std::abs(det);
  det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) +
          adztail * (bdx * cdy - cdx * bdy)) +
         (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) +
          bdztail * (cdx * ady - adx * cdy)) +
         (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) +
          cdztail * (adx * bdy - bdx * ady));
  if (det >= errbound or -det >= errbound)
    return det > 0.0 ? 1 : det < 0.0 ? -1 : 0;

  return Orient3DExact(a, b, c, d);
}

} // namespace

int
Orient3D(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d)
{
  const double adx = a.x - d.x;
  const double bdx = b.x - d.x;
  const double cdx = c.x - d.x;
  const double ady = a.y - d.y;
  const double bdy = b.y - d.y;
  const double cdy = c.y - d.y;
  const double adz = a.z - d.z;
  const double bdz = b.z - d.z;
  const double cdz = c.z - d.z;

  const double bdxcdy = bdx * cdy;
  const double cdxbdy = cdx * bdy;
  const double cdxady = cdx * ady;
  const double adxcdy = adx * cdy;
  const double adxbdy = adx * bdy;
  const double bdxady = bdx * ady;

  const double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
  const double permanent = (std::abs(bdxcdy) + std::abs(cdxbdy)) * std::abs(adz) +
                           (std::abs(cdxady) + std::abs(adxcdy)) * std::abs(bdz) +
                           (std::abs(adxbdy) + std::abs(bdxady)) * std::abs(cdz);
  // Shewchuk's bound A on the round-off of the floating-point evaluation.
  constexpr double epsilon = std::numeric_limits<double>::epsilon() / 2.0;
  constexpr double errbound = (7.0 + 56.0 * epsilon) * epsilon;
  if (det > errbound * permanent)
    return 1;
  if (-det > errbound * permanent)
    return -1;
  return Orient3DAdaptive(a, b, c, d, permanent);
}

} // namespace opensn
