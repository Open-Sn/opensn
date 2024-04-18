// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/acceleration.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{
namespace lbs
{

TwoGridCollapsedInfo
MakeTwoGridCollapsedInfo(const MultiGroupXS& xs, EnergyCollapseScheme scheme)
{
  const std::string fname = "lbs::acceleration::MakeTwoGridCollapsedInfo";

  const size_t num_groups = xs.NumGroups();
  const auto& sigma_t = xs.SigmaTotal();
  const auto& diffusion_coeff = xs.DiffusionCoefficient();

  // Make a Dense matrix from sparse transfer matrix
  if (xs.TransferMatrices().empty())
    throw std::logic_error(fname + ": list of scattering matrices empty.");

  const auto& isotropic_transfer_matrix = xs.TransferMatrix(0);

  MatDbl S(num_groups, VecDbl(num_groups, 0.0));
  for (int g = 0; g < num_groups; g++)
    for (const auto& [row_g, gprime, sigma] : isotropic_transfer_matrix.Row(g))
      S[g][gprime] = sigma;

  // Compiling the A and B matrices for different methods
  MatDbl A(num_groups, VecDbl(num_groups, 0.0));
  MatDbl B(num_groups, VecDbl(num_groups, 0.0));
  for (int g = 0; g < num_groups; g++)
  {
    if (scheme == EnergyCollapseScheme::JFULL)
    {
      A[g][g] = sigma_t[g] - S[g][g];
      for (int gp = 0; gp < g; gp++)
        B[g][gp] = S[g][gp];

      for (int gp = g + 1; gp < num_groups; gp++)
        B[g][gp] = S[g][gp];
    }
    else if (scheme == EnergyCollapseScheme::JPARTIAL)
    {
      A[g][g] = sigma_t[g];
      for (int gp = 0; gp < num_groups; gp++)
        B[g][gp] = S[g][gp];
    }
  } // for g

  // Correction for zero xs groups
  // Some cross-sections developed from monte-carlo
  // methods can result in some of the groups
  // having zero cross-sections. In that case
  // it will screw up the power iteration
  // initial guess of 1.0. Here we reset them
  for (int g = 0; g < num_groups; g++)
    if (sigma_t[g] < 1.0e-16)
      A[g][g] = 1.0;

  MatDbl Ainv = Inverse(A);
  MatDbl C = MatMul(Ainv, B);
  VecDbl E(num_groups, 1.0);

  double collapsed_D = 0.0;
  double collapsed_sig_a = 0.0;
  std::vector<double> spectrum(num_groups, 1.0);

  // Perform power iteration
  double rho = PowerIteration(C, E, 1000, 1.0e-12);

  // Compute two-grid diffusion quantities
  double sum = 0.0;
  for (int g = 0; g < num_groups; g++)
    sum += std::fabs(E[g]);

  for (int g = 0; g < num_groups; g++)
    spectrum[g] = std::fabs(E[g]) / sum;

  for (int g = 0; g < num_groups; ++g)
  {
    collapsed_D += diffusion_coeff[g] * spectrum[g];

    collapsed_sig_a += sigma_t[g] * spectrum[g];

    for (int gp = 0; gp < num_groups; ++gp)
      collapsed_sig_a -= S[g][gp] * spectrum[gp];
  }

  // Verbose output the spectrum
  log.Log0Verbose1() << "Fundamental eigen-value: " << rho;
  std::stringstream outstr;
  for (auto& xi : spectrum)
    outstr << xi << '\n';
  log.Log0Verbose1() << outstr.str();

  return {collapsed_D, collapsed_sig_a, spectrum};
}

} // namespace lbs
} // namespace opensn
