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
#include <string>
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

namespace
{

// The data that SLDFEsqQuadrature::LocallyRefine rebuilds: the directions, the harmonic selection,
// and the operators.
struct QuadratureState
{
  std::vector<double> phi;
  std::vector<double> theta;
  std::vector<double> weights;
  std::vector<double> omegas;
  std::vector<std::pair<unsigned int, int>> harmonics;
  unsigned int scattering_order = 0;
  std::vector<std::size_t> d2m_shape;
  std::vector<double> d2m;
  std::vector<std::size_t> m2d_shape;
  std::vector<double> m2d;
};

QuadratureState
StateOf(const AngularQuadrature& quadrature)
{
  QuadratureState state;
  for (const auto& abscissa : quadrature.GetAbscissae())
  {
    state.phi.push_back(abscissa.phi);
    state.theta.push_back(abscissa.theta);
  }
  state.weights = quadrature.GetWeights();
  for (const auto& omega : quadrature.GetOmegas())
    state.omegas.insert(state.omegas.end(), {omega.x, omega.y, omega.z});
  for (const auto& harmonic : quadrature.GetMomentToHarmonicsIndexMap())
    state.harmonics.emplace_back(harmonic.ell, harmonic.m);
  state.scattering_order = quadrature.GetScatteringOrder();
  const auto& d2m = quadrature.GetDiscreteToMomentOperator();
  state.d2m_shape = d2m.dimension();
  state.d2m.assign(d2m.begin(), d2m.end());
  const auto& m2d = quadrature.GetMomentToDiscreteOperator();
  state.m2d_shape = m2d.dimension();
  state.m2d.assign(m2d.begin(), m2d.end());
  return state;
}

// Compares two states exactly and names every part that differs.
::testing::AssertionResult
SameState(const QuadratureState& actual, const QuadratureState& expected)
{
  std::vector<std::string> changed;
  if (actual.phi != expected.phi or actual.theta != expected.theta)
    changed.emplace_back("abscissae");
  if (actual.weights != expected.weights)
    changed.emplace_back("weights");
  if (actual.omegas != expected.omegas)
    changed.emplace_back("directions");
  if (actual.harmonics != expected.harmonics)
    changed.emplace_back("harmonics");
  if (actual.scattering_order != expected.scattering_order)
    changed.emplace_back("scattering order");
  if (actual.d2m_shape != expected.d2m_shape or actual.d2m != expected.d2m)
    changed.emplace_back("discrete-to-moment operator");
  if (actual.m2d_shape != expected.m2d_shape or actual.m2d != expected.m2d)
    changed.emplace_back("moment-to-discrete operator");
  if (changed.empty())
    return ::testing::AssertionSuccess();

  auto failure = ::testing::AssertionFailure() << "changed:";
  for (std::size_t i = 0; i < changed.size(); ++i)
    failure << (i == 0 ? " " : ", ") << changed[i];
  failure << " (" << actual.omegas.size() / 3 << " directions and " << actual.harmonics.size()
          << " moments, expected " << expected.omegas.size() / 3 << " and "
          << expected.harmonics.size() << ")";
  return failure;
}

// A cone size of 0 splits no square, so the refinement rebuilds the set from its squares as they
// are. Returns how the refinement ended.
std::string
RefineWithoutSplitting(SLDFEsqQuadrature& quadrature)
{
  try
  {
    quadrature.LocallyRefine(Vector3(0.0, 0.0, 1.0), 0.0);
  }
  catch (const std::invalid_argument&)
  {
    return "rejected";
  }
  catch (const std::exception& e)
  {
    return std::string("failed (") + e.what() + ")";
  }
  return "accepted";
}

// The squares of an SLDFE set are not public; only a later refinement reads them. If they are
// unchanged, a refinement that splits none of them ends as on a new set and leaves the same state.
::testing::AssertionResult
RefinesLikeNewSet(SLDFEsqQuadrature& quadrature, SLDFEsqQuadrature& new_set)
{
  const auto outcome = RefineWithoutSplitting(quadrature);
  const auto expected = RefineWithoutSplitting(new_set);
  if (outcome != expected)
    return ::testing::AssertionFailure()
           << "a refinement that splits no square was " << outcome << ", on a new set " << expected;
  return SameState(StateOf(quadrature), StateOf(new_set));
}

// GALERKIN_ONE builds the discrete-to-moment operator from the inverse of the moment-to-discrete
// operator, so for a set with n directions both are n x n and sum_d D2M(d, m) M2D(d, k) = delta_mk.
::testing::AssertionResult
HasInverseOperators(const AngularQuadrature& quadrature, std::size_t num_directions)
{
  if (quadrature.GetNumAngles() != num_directions or quadrature.GetNumMoments() != num_directions)
    return ::testing::AssertionFailure()
           << quadrature.GetNumAngles() << " directions and " << quadrature.GetNumMoments()
           << " moments, expected " << num_directions << " of each";

  const auto& d2m = quadrature.GetDiscreteToMomentOperator();
  const auto& m2d = quadrature.GetMomentToDiscreteOperator();
  const std::vector<std::size_t> shape = {num_directions, num_directions};
  if (d2m.dimension() != shape or m2d.dimension() != shape)
    return ::testing::AssertionFailure()
           << "operator shapes (" << d2m.dimension()[0] << ", " << d2m.dimension()[1] << ") and ("
           << m2d.dimension()[0] << ", " << m2d.dimension()[1] << "), expected (" << num_directions
           << ", " << num_directions << ")";

  double max_error = 0.0;
  for (std::size_t m = 0; m < num_directions; ++m)
    for (std::size_t k = 0; k < num_directions; ++k)
    {
      double sum = 0.0;
      for (std::size_t d = 0; d < num_directions; ++d)
        sum += d2m(d, m) * m2d(d, k);
      max_error = std::max(max_error, std::fabs(sum - (m == k ? 1.0 : 0.0)));
    }
  if (max_error > 1.0e-9)
    return ::testing::AssertionFailure() << "max |D2M^T M2D - I| = " << max_error;

  return ::testing::AssertionSuccess();
}

} // namespace

TEST(GalerkinOneQuadratureTest, LocallyRefinedSLDFEThrows)
{
  // A locally refined level-0 set has a direction count that no rule describes. The rejection
  // must leave the set as it was.
  const auto method = OperatorConstructionMethod::GALERKIN_ONE;
  SLDFEsqQuadrature3DXYZ sldfe(0, 0, false, method);
  const auto before = StateOf(sldfe);
  EXPECT_TRUE(IsRejected([&sldfe] { sldfe.LocallyRefine(Vector3(0.0, 0.0, 1.0), 1.0); }));
  EXPECT_TRUE(SameState(StateOf(sldfe), before));
  SLDFEsqQuadrature3DXYZ new_sldfe(0, 0, false, method);
  EXPECT_TRUE(RefinesLikeNewSet(sldfe, new_sldfe));

  // The same in 2D, where the set keeps only the directions with a positive z component.
  SLDFEsqQuadrature2DXY sldfe_2d(0, 0, false, method);
  const auto before_2d = StateOf(sldfe_2d);
  EXPECT_TRUE(IsRejected([&sldfe_2d] { sldfe_2d.LocallyRefine(Vector3(0.0, 0.0, 1.0), 1.0); }));
  EXPECT_TRUE(SameState(StateOf(sldfe_2d), before_2d));
  SLDFEsqQuadrature2DXY new_sldfe_2d(0, 0, false, method);
  EXPECT_TRUE(RefinesLikeNewSet(sldfe_2d, new_sldfe_2d));
}

TEST(GalerkinOneQuadratureTest, LocallyRefinedSLDFEWithRuleBuilds)
{
  // A cone size of pi splits every square of a level-0 set, which then has the direction count of
  // level 1 (384 in 3D, 192 in 2D). A level-0 set takes its rule from its direction count, so the
  // level-1 rule applies, and both operators must be rebuilt for the refined directions.
  const auto method = OperatorConstructionMethod::GALERKIN_ONE;
  SLDFEsqQuadrature3DXYZ sldfe(0, 0, false, method);
  EXPECT_NO_THROW(sldfe.LocallyRefine(Vector3(0.0, 0.0, 1.0), M_PI));
  EXPECT_TRUE(HasInverseOperators(sldfe, 384));

  SLDFEsqQuadrature2DXY sldfe_2d(0, 0, false, method);
  EXPECT_NO_THROW(sldfe_2d.LocallyRefine(Vector3(0.0, 0.0, 1.0), M_PI));
  EXPECT_TRUE(HasInverseOperators(sldfe_2d, 192));
}

TEST(GalerkinOneQuadratureTest, LocallyRefinedSLDFEWithMismatchedRuleThrows)
{
  // Above level 0 an SLDFE set takes its rule from its level, so a refined level-1 set still
  // selects with the level-1 rule, which gives fewer harmonics than the set has directions (384
  // for 672 in 3D, 192 for 480 in 2D). The rejection must leave the set as it was.
  const auto method = OperatorConstructionMethod::GALERKIN_ONE;
  SLDFEsqQuadrature3DXYZ sldfe(1, 0, false, method);
  const auto before = StateOf(sldfe);
  EXPECT_TRUE(IsRejected([&sldfe] { sldfe.LocallyRefine(Vector3(0.0, 0.0, 1.0), 1.0); }));
  EXPECT_TRUE(SameState(StateOf(sldfe), before));

  SLDFEsqQuadrature2DXY sldfe_2d(1, 0, false, method);
  const auto before_2d = StateOf(sldfe_2d);
  EXPECT_TRUE(IsRejected([&sldfe_2d] { sldfe_2d.LocallyRefine(Vector3(0.0, 0.0, 1.0), 1.0); }));
  EXPECT_TRUE(SameState(StateOf(sldfe_2d), before_2d));
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
