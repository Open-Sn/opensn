// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/lebedev_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/data_types/vector3.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace opensn;

// GALERKIN_ONE needs a square moment-to-discrete operator, so it selects as many harmonics as
// there are directions. For Lebedev and SLDFE sets that selection comes from tabulated rules.
// A set without a tabulated rule cannot be built with GALERKIN_ONE and must be rejected with
// std::invalid_argument. If it is not rejected, the harmonic selection does not terminate and
// allocates without bound, so run these tests under a memory limit (for example
// `ulimit -v 4194304`). Any other exception, such as the std::bad_alloc that ends such a run,
// counts as a failure, not as a rejection.

namespace
{

template <typename Construct>
::testing::AssertionResult
IsRejected(Construct construct)
{
  try
  {
    construct();
  }
  catch (const std::invalid_argument&)
  {
    return ::testing::AssertionSuccess();
  }
  catch (const std::bad_alloc&)
  {
    return ::testing::AssertionFailure()
           << "std::bad_alloc: the harmonic selection ran out of memory instead of rejecting "
              "the set";
  }
  catch (const std::exception& e)
  {
    return ::testing::AssertionFailure()
           << "the set was rejected with an exception other than std::invalid_argument: "
           << e.what();
  }
  catch (...)
  {
    return ::testing::AssertionFailure() << "the set was rejected with a non-standard exception";
  }
  return ::testing::AssertionFailure() << "the set was accepted";
}

} // namespace

TEST(GalerkinOneQuadratureTest, TabulatedRulesGiveSquareOperators)
{
  const LebedevQuadrature3DXYZ lebedev(3, 0, false, OperatorConstructionMethod::GALERKIN_ONE);
  EXPECT_EQ(lebedev.GetNumMoments(), lebedev.GetNumAngles());

  const SLDFEsqQuadrature3DXYZ sldfe(0, 0, false, OperatorConstructionMethod::GALERKIN_ONE);
  EXPECT_EQ(sldfe.GetNumMoments(), sldfe.GetNumAngles());
}

TEST(GalerkinOneQuadratureTest, LebedevOrderWithoutRuleThrows)
{
  // Lebedev points exist up to order 131; the harmonic rules stop at order 101.
  EXPECT_TRUE(IsRejected(
    [] { LebedevQuadrature3DXYZ(107, 0, false, OperatorConstructionMethod::GALERKIN_ONE); }));
  EXPECT_TRUE(IsRejected(
    [] { LebedevQuadrature2DXY(107, 0, false, OperatorConstructionMethod::GALERKIN_ONE); }));
}

TEST(GalerkinOneQuadratureTest, SLDFELevelWithoutRuleThrows)
{
  // The SLDFE harmonic rules cover refinement levels 0 to 3.
  EXPECT_TRUE(IsRejected(
    [] { SLDFEsqQuadrature3DXYZ(4, 0, false, OperatorConstructionMethod::GALERKIN_ONE); }));
  EXPECT_TRUE(IsRejected(
    [] { SLDFEsqQuadrature2DXY(4, 0, false, OperatorConstructionMethod::GALERKIN_ONE); }));
}

TEST(GalerkinOneQuadratureTest, LocallyRefinedSLDFEThrows)
{
  // A locally refined level-0 set has a direction count that no rule describes.
  SLDFEsqQuadrature3DXYZ sldfe(0, 0, false, OperatorConstructionMethod::GALERKIN_ONE);
  EXPECT_TRUE(IsRejected([&sldfe] { sldfe.LocallyRefine(Vector3(0.0, 0.0, 1.0), 1.0); }));
}

namespace
{

// The standard harmonic selection for a scattering order: in 3D every (ell, m) with |m| <= ell,
// in 2D only m = -ell, -ell + 2, ..., ell (the harmonics that are even in the polar cosine).
std::vector<std::pair<unsigned int, int>>
StandardHarmonics(unsigned int dimension, unsigned int scattering_order)
{
  std::vector<std::pair<unsigned int, int>> harmonics;
  const int m_step = (dimension == 2) ? 2 : 1;
  const int max_ell = static_cast<int>(scattering_order);
  for (int ell = 0; ell <= max_ell; ++ell)
    for (int m = -ell; m <= ell; m += m_step)
      harmonics.emplace_back(static_cast<unsigned int>(ell), m);
  return harmonics;
}

// GALERKIN_THREE orthogonalizes the harmonics under the quadrature weights and scales them to the
// norm 1/(2 ell + 1), so sum_n D2M(n, m) M2D(n, k) = delta_mk / W, with W the weight sum.
::testing::AssertionResult
BuildsWithStandardSelection(const AngularQuadrature& quadrature,
                            unsigned int dimension,
                            unsigned int scattering_order)
{
  if (quadrature.GetScatteringOrder() != scattering_order)
    return ::testing::AssertionFailure() << "scattering order " << quadrature.GetScatteringOrder()
                                         << ", expected " << scattering_order;

  const auto expected = StandardHarmonics(dimension, scattering_order);
  const auto& selected = quadrature.GetMomentToHarmonicsIndexMap();
  if (selected.size() != expected.size())
    return ::testing::AssertionFailure()
           << selected.size() << " moments, expected " << expected.size();
  for (std::size_t i = 0; i < expected.size(); ++i)
    if (selected[i].ell != expected[i].first or selected[i].m != expected[i].second)
      return ::testing::AssertionFailure()
             << "moment " << i << " is (" << selected[i].ell << ", " << selected[i].m
             << "), expected (" << expected[i].first << ", " << expected[i].second << ")";

  const std::size_t num_angles = quadrature.GetNumAngles();
  const std::size_t num_moments = expected.size();
  const auto& d2m = quadrature.GetDiscreteToMomentOperator();
  const auto& m2d = quadrature.GetMomentToDiscreteOperator();
  const std::vector<std::size_t> shape = {num_angles, num_moments};
  if (d2m.dimension() != shape or m2d.dimension() != shape)
    return ::testing::AssertionFailure()
           << "operator shapes differ from (" << num_angles << ", " << num_moments << ")";

  const double inverse_weight_sum = 1.0 / quadrature.GetWeightSum();
  double max_error = 0.0;
  for (std::size_t m = 0; m < num_moments; ++m)
    for (std::size_t k = 0; k < num_moments; ++k)
    {
      double sum = 0.0;
      for (std::size_t n = 0; n < num_angles; ++n)
        sum += d2m(n, m) * m2d(n, k);
      max_error = std::max(max_error, std::fabs(sum - (m == k ? inverse_weight_sum : 0.0)));
    }
  if (max_error > 1.0e-12)
    return ::testing::AssertionFailure() << "max |D2M^T M2D - I / W| = " << max_error;

  return ::testing::AssertionSuccess();
}

template <typename Construct>
::testing::AssertionResult
BuildsWithGalerkinThree(Construct construct, unsigned int dimension, unsigned int scattering_order)
{
  std::unique_ptr<AngularQuadrature> quadrature;
  try
  {
    quadrature = construct();
  }
  catch (const std::exception& e)
  {
    return ::testing::AssertionFailure() << "the set was rejected: " << e.what();
  }
  return BuildsWithStandardSelection(*quadrature, dimension, scattering_order);
}

} // namespace

TEST(GalerkinOneQuadratureTest, SetsWithoutRuleBuildWithGalerkinThree)
{
  // Without a rule GALERKIN_THREE falls back to the standard selection for the given scattering
  // order, so the sets that GALERKIN_ONE must reject still build with GALERKIN_THREE.
  constexpr unsigned int L = 3;
  const auto method = OperatorConstructionMethod::GALERKIN_THREE;
  EXPECT_TRUE(BuildsWithGalerkinThree(
    [&] { return std::make_unique<LebedevQuadrature3DXYZ>(107, L, false, method); }, 3, L));
  EXPECT_TRUE(BuildsWithGalerkinThree(
    [&] { return std::make_unique<LebedevQuadrature2DXY>(107, L, false, method); }, 2, L));
  EXPECT_TRUE(BuildsWithGalerkinThree(
    [&] { return std::make_unique<SLDFEsqQuadrature3DXYZ>(4, L, false, method); }, 3, L));
  EXPECT_TRUE(BuildsWithGalerkinThree(
    [&] { return std::make_unique<SLDFEsqQuadrature2DXY>(4, L, false, method); }, 2, L));
  EXPECT_TRUE(BuildsWithGalerkinThree(
    [&]
    {
      auto sldfe = std::make_unique<SLDFEsqQuadrature3DXYZ>(0, L, false, method);
      sldfe->LocallyRefine(Vector3(0.0, 0.0, 1.0), 1.0);
      return sldfe;
    },
    3,
    L));
}
