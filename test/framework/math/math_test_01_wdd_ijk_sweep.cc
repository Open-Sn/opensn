#include <array>
#include <cstddef>
#include <cmath>

#include "framework/data_types/ndarray.h"
#include "framework/math/quadratures/angular_product_quadrature.h"
#include "framework/math/math_range.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "framework/console/console.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock chi_math_Test01_WDD_IJK_Sweep(const InputParameters& params);

RegisterWrapperFunction(chi_unit_tests,
                        chi_math_Test01_WDD_IJK_Sweep,
                        nullptr,
                        chi_math_Test01_WDD_IJK_Sweep);

typedef NDArray<double> IJKArrayDbl;

IJKArrayDbl
WDD_IJK_Sweep2(const std::array<size_t, 3>& mesh_divs,
               const std::array<double, 3>& mesh_lengths,
               const std::array<double, 6>& bcs,
               const IJKArrayDbl& sigma_t,
               const IJKArrayDbl& q,
               const AngularQuadrature& quad,
               bool verbose = false)
{
  const int Nx = static_cast<int>(mesh_divs[0]);
  const int Ny = static_cast<int>(mesh_divs[1]);
  const int Nz = static_cast<int>(mesh_divs[2]);

  const double dx = mesh_lengths[0] / Nx;
  const double dy = mesh_lengths[1] / Ny;
  const double dz = mesh_lengths[2] / Nz;

  IJKArrayDbl phi_0(mesh_divs);
  phi_0.set(0.0);

  auto iorder = Range<int>(0, Nx);
  auto jorder = Range<int>(0, Ny);
  auto korder = Range<int>(0, Nz);

  const auto& D2M = quad.GetDiscreteToMomentOperator();

  int n = 0;
  for (const auto& omega_n : quad.omegas_)
  {
    if (opensn::mpi.location_id == 0 and verbose)
      std::cout << "Sweep angle " << n << " " << omega_n.PrintStr() << std::endl;

    // Determine sweep ordering
    if (omega_n.x > 0.0) iorder = Range<int>(0, Nx);
    else
      iorder = Range<int>(Nx - 1, -1, -1);

    if (omega_n.y > 0.0) jorder = Range<int>(0, Ny);
    else
      jorder = Range<int>(Ny - 1, -1, -1);

    if (omega_n.z > 0.0) korder = Range<int>(0, Nz);
    else
      korder = Range<int>(Nz - 1, -1, -1);

    // Sweep cells
    IJKArrayDbl psi_ds_x(mesh_divs);
    IJKArrayDbl psi_ds_y(mesh_divs);
    IJKArrayDbl psi_ds_z(mesh_divs);
    for (auto k : korder)
      for (auto j : jorder)
        for (auto i : iorder)
        {
          double psi_us_x = (omega_n.x > 0.0) ? bcs[0] : bcs[1];
          double psi_us_y = (omega_n.y > 0.0) ? bcs[2] : bcs[3];
          double psi_us_z = (omega_n.z > 0.0) ? bcs[4] : bcs[5];

          if (omega_n.x > 0.0 and i > 0) psi_us_x = psi_ds_x(i - 1, j, k);
          if (omega_n.x < 0.0 and i < (Nx - 1)) psi_us_x = psi_ds_x(i + 1, j, k);

          if (omega_n.y > 0.0 and j > 0) psi_us_y = psi_ds_y(i, j - 1, k);
          if (omega_n.y < 0.0 and j < (Ny - 1)) psi_us_y = psi_ds_y(i, j + 1, k);

          if (omega_n.z > 0.0 and k > 0) psi_us_z = psi_ds_z(i, j, k - 1);
          if (omega_n.z < 0.0 and k < (Nz - 1)) psi_us_z = psi_ds_z(i, j, k + 1);

          double rhs = q(i, j, k) / (4.0 * M_PI);
          if (Nx > 1) rhs += 2.0 * std::fabs(omega_n.x) * psi_us_x / dx;
          if (Ny > 1) rhs += 2.0 * std::fabs(omega_n.y) * psi_us_y / dy;
          if (Nz > 1) rhs += 2.0 * std::fabs(omega_n.z) * psi_us_z / dz;

          double lhs = sigma_t(i, j, k);
          if (Nx > 1) lhs += 2.0 * std::fabs(omega_n.x) / dx;
          if (Ny > 1) lhs += 2.0 * std::fabs(omega_n.y) / dy;
          if (Nz > 1) lhs += 2.0 * std::fabs(omega_n.z) / dz;

          double psi_ijk = rhs / lhs;

          phi_0(i, j, k) += D2M[0][n] * psi_ijk;

          psi_ds_x(i, j, k) = 2.0 * psi_ijk - psi_us_x;
          psi_ds_y(i, j, k) = 2.0 * psi_ijk - psi_us_y;
          psi_ds_z(i, j, k) = 2.0 * psi_ijk - psi_us_z;
        }
    ++n;
  } // for omega

  return phi_0;
}

ParameterBlock
chi_math_Test01_WDD_IJK_Sweep(const InputParameters&)
{
  opensn::log.Log() << "GOLD_BEGIN";
  bool verbose = true;
  const std::array<size_t, 3> mesh_divisions = {1, 1, 10};
  const std::array<double, 3> mesh_lengths = {1.0, 1.0, 10.0};
  const std::array<double, 6> bcs = {0.0, 0.0, 0.0, 0.0, 0.5, 0.0};

  IJKArrayDbl sigma_t(mesh_divisions, 0.2);
  IJKArrayDbl q(mesh_divisions, 0.0);

  //  sigma_t.Set(0.2);
  //  q.Set(0.0);

  auto pquad = std::make_shared<AngularQuadratureProdGL>(1, verbose);

  pquad->BuildDiscreteToMomentOperator(0, 1);

  auto phi = WDD_IJK_Sweep2(mesh_divisions, mesh_lengths, bcs, sigma_t, q, *pquad, verbose);

  if (opensn::mpi.location_id == 0 and verbose)
  {
    std::cout << "order:\n";
    for (auto i : phi)
      std::cout << i << "\n";
  }

  opensn::log.Log() << "GOLD_END";
  return ParameterBlock();
}

} //  namespace unit_tests
