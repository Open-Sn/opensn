// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "gmock/gmock.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/data_types/ndarray.h"
#include "framework/data_types/range.h"
#include <array>
#include <cstddef>
#include <cmath>

using namespace opensn;

namespace
{

NDArray<double, 3>
WDD_IJK_Sweep2(const std::array<size_t, 3>& mesh_divs,
               const std::array<double, 3>& mesh_lengths,
               const std::array<double, 6>& bcs,
               const NDArray<double, 3>& sigma_t,
               const NDArray<double, 3>& q,
               const AngularQuadrature& quad,
               bool verbose = false)
{
  const int Nx = static_cast<int>(mesh_divs[0]);
  const int Ny = static_cast<int>(mesh_divs[1]);
  const int Nz = static_cast<int>(mesh_divs[2]);

  const double dx = mesh_lengths[0] / Nx;
  const double dy = mesh_lengths[1] / Ny;
  const double dz = mesh_lengths[2] / Nz;

  NDArray<double, 3> phi_0(mesh_divs, 0.);

  auto iorder = Range<int>(0, Nx);
  auto jorder = Range<int>(0, Ny);
  auto korder = Range<int>(0, Nz);

  const auto& D2M = quad.GetDiscreteToMomentOperator();

  int n = 0;
  for (const auto& omega_n : quad.omegas)
  {
    // Determine sweep ordering
    if (omega_n.x > 0.0)
      iorder = Range<int>(0, Nx);
    else
      iorder = Range<int>(Nx - 1, -1, -1);

    if (omega_n.y > 0.0)
      jorder = Range<int>(0, Ny);
    else
      jorder = Range<int>(Ny - 1, -1, -1);

    if (omega_n.z > 0.0)
      korder = Range<int>(0, Nz);
    else
      korder = Range<int>(Nz - 1, -1, -1);

    // Sweep cells
    NDArray<double, 3> psi_ds_x(mesh_divs);
    NDArray<double, 3> psi_ds_y(mesh_divs);
    NDArray<double, 3> psi_ds_z(mesh_divs);
    for (auto k : korder)
      for (auto j : jorder)
        for (auto i : iorder)
        {
          double psi_us_x = (omega_n.x > 0.0) ? bcs[0] : bcs[1];
          double psi_us_y = (omega_n.y > 0.0) ? bcs[2] : bcs[3];
          double psi_us_z = (omega_n.z > 0.0) ? bcs[4] : bcs[5];

          if (omega_n.x > 0.0 and i > 0)
            psi_us_x = psi_ds_x(i - 1, j, k);
          if (omega_n.x < 0.0 and i < (Nx - 1))
            psi_us_x = psi_ds_x(i + 1, j, k);

          if (omega_n.y > 0.0 and j > 0)
            psi_us_y = psi_ds_y(i, j - 1, k);
          if (omega_n.y < 0.0 and j < (Ny - 1))
            psi_us_y = psi_ds_y(i, j + 1, k);

          if (omega_n.z > 0.0 and k > 0)
            psi_us_z = psi_ds_z(i, j, k - 1);
          if (omega_n.z < 0.0 and k < (Nz - 1))
            psi_us_z = psi_ds_z(i, j, k + 1);

          double rhs = q(i, j, k) / (4.0 * M_PI);
          if (Nx > 1)
            rhs += 2.0 * std::fabs(omega_n.x) * psi_us_x / dx;
          if (Ny > 1)
            rhs += 2.0 * std::fabs(omega_n.y) * psi_us_y / dy;
          if (Nz > 1)
            rhs += 2.0 * std::fabs(omega_n.z) * psi_us_z / dz;

          double lhs = sigma_t(i, j, k);
          if (Nx > 1)
            lhs += 2.0 * std::fabs(omega_n.x) / dx;
          if (Ny > 1)
            lhs += 2.0 * std::fabs(omega_n.y) / dy;
          if (Nz > 1)
            lhs += 2.0 * std::fabs(omega_n.z) / dz;

          double psi_ijk = rhs / lhs;

          phi_0(i, j, k) += D2M(0, n) * psi_ijk;

          psi_ds_x(i, j, k) = 2.0 * psi_ijk - psi_us_x;
          psi_ds_y(i, j, k) = 2.0 * psi_ijk - psi_us_y;
          psi_ds_z(i, j, k) = 2.0 * psi_ijk - psi_us_z;
        }
    ++n;
  } // for omega

  return phi_0;
}

} // namespace

TEST(MathTest, WDD_IJK_Sweep)
{
  bool verbose = false;
  const std::array<size_t, 3> mesh_divisions = {1, 1, 10};
  const std::array<double, 3> mesh_lengths = {1.0, 1.0, 10.0};
  const std::array<double, 6> bcs = {0.0, 0.0, 0.0, 0.0, 0.5, 0.0};

  NDArray<double, 3> sigma_t(mesh_divisions, 0.2);
  NDArray<double, 3> q(mesh_divisions, 0.0);

  auto pquad = std::make_shared<GLProductQuadrature1DSlab>(2, 0, verbose);
  auto phi = WDD_IJK_Sweep2(mesh_divisions, mesh_lengths, bcs, sigma_t, q, *pquad, verbose);

  {
    ASSERT_EQ(pquad->omegas.size(), 2);
    EXPECT_NEAR(pquad->omegas[0].x, 0.816497, 1e-5);
    EXPECT_NEAR(pquad->omegas[0].y, 0, 1e-5);
    EXPECT_NEAR(pquad->omegas[0].z, 0.57735, 1e-5);
    EXPECT_NEAR(pquad->omegas[1].x, 0.816497, 1e-5);
    EXPECT_NEAR(pquad->omegas[1].y, 0., 1e-5);
    EXPECT_NEAR(pquad->omegas[1].z, -0.57735, 1e-5);
  }

  {
    auto* vals = phi.data();
    ASSERT_EQ(phi.size(), 10);
    EXPECT_NEAR(vals[0], 0.21309152596827, 1e-6);
    EXPECT_NEAR(vals[1], 0.15017266359160633, 1e-6);
    EXPECT_NEAR(vals[2], 0.10583078021596619, 1e-6);
    EXPECT_NEAR(vals[3], 0.074582374805187826, 1e-6);
    EXPECT_NEAR(vals[4], 0.052560252875463634, 1e-6);
    EXPECT_NEAR(vals[5], 0.037041270476724626, 1e-6);
    EXPECT_NEAR(vals[6], 0.026104167930776304, 1e-6);
    EXPECT_NEAR(vals[7], 0.018396388221625048, 1e-6);
    EXPECT_NEAR(vals[8], 0.012964563756386171, 1e-6);
    EXPECT_NEAR(vals[9], 0.0091364791372981619, 1e-6);
  }
}
