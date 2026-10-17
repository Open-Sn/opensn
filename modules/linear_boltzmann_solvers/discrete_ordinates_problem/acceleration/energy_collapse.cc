// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/acceleration.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <cmath>

namespace opensn
{

TwoGridCollapsedInfo
MakeTwoGridCollapsedInfo(const MultiGroupXS& xs,
                         EnergyCollapseScheme scheme,
                         const unsigned int first_group,
                         const unsigned int last_group,
                         const double time_absorption_scale)
{
  const std::string fname = "acceleration::MakeTwoGridCollapsedInfo";

  const auto num_groups = xs.GetNumGroups();
  if (first_group > last_group or last_group >= num_groups)
    throw std::logic_error(fname + ": invalid group range [" + std::to_string(first_group) + ", " +
                           std::to_string(last_group) + "] for a cross section with " +
                           std::to_string(num_groups) + " groups.");
  // The two-grid error mode lives in the groups of the accelerated groupset. Groups outside it are
  // not iterated with the groupset, so they are excluded from the collapse; scattering from the
  // groupset to them is removal (it is part of the total cross section).
  const unsigned int n = last_group - first_group + 1;
  const auto& diffusion_coeff_xs = xs.GetDiffusionCoefficient();

  // Include the time absorption 1/(v theta dt) of the theta scheme in the total cross section
  // and the diffusion coefficient.
  std::vector<double> sigma_t(n);
  std::vector<double> diffusion_coeff(n);
  for (unsigned int i = 0; i < n; ++i)
  {
    sigma_t[i] = xs.GetSigmaTotal()[first_group + i];
    diffusion_coeff[i] = diffusion_coeff_xs[first_group + i];
  }
  if (time_absorption_scale > 0.0)
  {
    const auto& inv_velocity = xs.GetInverseVelocity();
    if (inv_velocity.size() < num_groups)
      throw std::logic_error(fname + ": time-dependent two-grid acceleration requires inverse "
                                     "velocities.");
    for (unsigned int i = 0; i < n; ++i)
    {
      const double tau = inv_velocity[first_group + i] * time_absorption_scale;
      sigma_t[i] += tau;
      diffusion_coeff[i] =
        AddTimeAbsorptionToDiffusionCoefficient(diffusion_coeff_xs[first_group + i], tau);
    }
  }

  // Make a Dense matrix from sparse transfer matrix
  if (xs.GetTransferMatrices().empty())
    throw std::logic_error(fname + ": list of scattering matrices empty.");

  const auto& isotropic_transfer_matrix = xs.GetTransferMatrix(0);

  DenseMatrix<double> S(n, n, 0.0);
  for (unsigned int i = 0; i < n; ++i)
    for (const auto& [row_g, gprime, sigma] : isotropic_transfer_matrix.Row(first_group + i))
      if (gprime >= first_group and gprime <= last_group)
        S(i, gprime - first_group) = sigma;

  // Compiling the A and B matrices for different methods
  DenseMatrix<double> A(n, n, 0.0);
  DenseMatrix<double> B(n, n, 0.0);
  for (unsigned int g = 0; g < n; ++g)
  {
    if (scheme == EnergyCollapseScheme::JFULL)
    {
      A(g, g) = sigma_t[g] - S(g, g);
      for (unsigned int gp = 0; gp < g; ++gp)
        B(g, gp) = S(g, gp);

      for (unsigned int gp = g + 1; gp < n; ++gp)
        B(g, gp) = S(g, gp);
    }
    else if (scheme == EnergyCollapseScheme::JPARTIAL)
    {
      A(g, g) = sigma_t[g];
      for (unsigned int gp = 0; gp < n; ++gp)
        B(g, gp) = S(g, gp);
    }
  } // for g

  // Correction for zero xs groups
  // Some cross sections developed from monte-carlo
  // methods can result in some of the groups
  // having zero cross sections. In that case
  // it will screw up the power iteration
  // initial guess of 1.0. Here we reset them
  for (unsigned int g = 0; g < n; ++g)
    if (sigma_t[g] < 1.0e-16)
      A(g, g) = 1.0;

  auto Ainv = Inverse(A);
  auto C = Mult(Ainv, B);
  Vector<double> E(n, 1.0);

  double collapsed_D = 0.0;
  double collapsed_sig_a = 0.0;
  std::vector<double> block_spectrum(n, 1.0);

  // Perform power iteration. Without group-to-group scattering within the groupset there is no
  // two-grid error mode coupling the groups; a flat spectrum is used instead. The same holds when
  // the groupset has only downscatter and no lagged within-group scattering: C is then strictly
  // lower triangular (nilpotent), the iteration error vanishes after n iterations, and there is no
  // dominant eigenvector.
  bool has_group_coupling = false;
  bool has_upscatter_or_diagonal = false;
  for (unsigned int g = 0; g < n; ++g)
    for (unsigned int gp = 0; gp < n; ++gp)
    {
      if (g != gp)
        has_group_coupling = has_group_coupling or C(g, gp) != 0.0;
      if (gp >= g)
        has_upscatter_or_diagonal = has_upscatter_or_diagonal or C(g, gp) != 0.0;
    }
  const bool has_dominant_mode = has_group_coupling and has_upscatter_or_diagonal;

  // The spectral shape is the eigenvector of the largest real (Perron) eigenvalue rho of the
  // non-negative iteration matrix C. Plain power iteration finds the eigenvalue of largest
  // magnitude, which need not be unique: Jacobi iteration matrices can have -rho as an eigenvalue
  // as well (e.g. two groups coupled by down- and upscatter), and power iteration then oscillates
  // and returns a meaningless vector. Iterating with C + s I, s >= rho, makes rho + s strictly
  // dominant. The Gershgorin bound s = max_g sum_g' |C(g,g')| satisfies s >= rho.
  double rho = 0.0;
  double sum = 0.0;
  if (has_dominant_mode)
  {
    double shift = 0.0;
    for (unsigned int g = 0; g < n; ++g)
    {
      double row_sum = 0.0;
      for (unsigned int gp = 0; gp < n; ++gp)
        row_sum += std::fabs(C(g, gp));
      shift = std::max(shift, row_sum);
    }

    constexpr int max_iterations = 100000;
    constexpr double tolerance = 1.0e-13;
    for (unsigned int g = 0; g < n; ++g)
      E(g) = 1.0 / n;
    bool converged = false;
    for (int it = 0; it < max_iterations and not converged; ++it)
    {
      Vector<double> E_new = Mult(C, E);
      double norm = 0.0;
      for (unsigned int g = 0; g < n; ++g)
      {
        E_new(g) += shift * E(g);
        norm += std::fabs(E_new(g));
      }
      if (not(norm > 0.0) or not std::isfinite(norm))
        break;
      double change = 0.0;
      for (unsigned int g = 0; g < n; ++g)
      {
        E_new(g) /= norm;
        change = std::max(change, std::fabs(E_new(g) - E(g)));
      }
      rho = norm - shift; // E has unit L1 norm and is non-negative
      E = E_new;
      converged = change <= tolerance;
    }
    if (not converged)
      log.Log0Warning() << fname
                        << ": two-grid spectrum did not converge; the TGDSA correction may be "
                           "less effective.";
    for (unsigned int g = 0; g < n; ++g)
      sum += std::fabs(E(g));
  }

  // Compute two-grid diffusion quantities
  if (has_dominant_mode and std::isfinite(sum) and sum > 0.0)
    for (unsigned int g = 0; g < n; ++g)
      block_spectrum[g] = std::fabs(E(g)) / sum;
  else
  {
    log.Log0Warning() << fname
                      << ": material has no lagged group-to-group scattering mode in groups ["
                      << first_group << ", " << last_group
                      << "]; using a flat two-grid spectrum for it.";
    std::fill(block_spectrum.begin(), block_spectrum.end(), 1.0 / n);
  }

  for (unsigned int g = 0; g < n; ++g)
  {
    collapsed_D += diffusion_coeff[g] * block_spectrum[g];

    collapsed_sig_a += sigma_t[g] * block_spectrum[g];

    for (unsigned int gp = 0; gp < n; ++gp)
      collapsed_sig_a -= S(g, gp) * block_spectrum[gp];
  }

  // The spectrum is indexed by global group and is zero outside the groupset.
  std::vector<double> spectrum(num_groups, 0.0);
  for (unsigned int g = 0; g < n; ++g)
    spectrum[first_group + g] = block_spectrum[g];

  // Verbose output the spectrum
  log.Log0Verbose1() << "Fundamental eigen-value: " << rho;
  std::stringstream outstr;
  for (auto& xi : block_spectrum)
    outstr << xi << '\n';
  log.Log0Verbose1() << outstr.str();

  TwoGridCollapsedInfo tgci{};
  tgci.collapsed_D = collapsed_D;
  tgci.collapsed_sig_a = collapsed_sig_a;
  tgci.spectrum = spectrum;
  return tgci;
}

} // namespace opensn
