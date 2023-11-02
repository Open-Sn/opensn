#include "mg_diffusion_solver.h"
#include "chi_runtime.h"
#include "chi_log.h"
#include "physics/PhysicsMaterial/chi_physicsmaterial.h"

void mg_diffusion::Solver::Compute_TwoGrid_Params()
{
  // loop over all materials
  for (const auto &mat_id_xs: matid_to_xs_map) {

    // get the P0 transfer matrix and total XS
    const auto &isotropic_transfer_matrix = mat_id_xs.second->TransferMatrix(0);
    const auto &sigma_t = mat_id_xs.second->SigmaTotal();
    const auto &diffusion_coeff = mat_id_xs.second->DiffusionCoefficient();

    // put P0 transfer matrix in nicer form
    MatDbl S(num_groups_, VecDbl(num_groups_, 0.0));
    for (unsigned int g = 0; g < num_groups_; ++g)
      for (const auto &[row_g, gprime, sigma]: isotropic_transfer_matrix.Row(g))
        S[g][gprime] = sigma;

    // (L+D) e_new = -U e_old
    // original matrix = diag(total) - scattering
    // so L+D = diag(removal) - tril(scattering)
    // and U = -triu(scattering)
    MatDbl A(num_groups_, VecDbl(num_groups_, 0.0));
    MatDbl B(num_groups_, VecDbl(num_groups_, 0.0));
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      A[g][g] = sigma_t[g] - S[g][g];
      for (unsigned int gp = 0; gp < g; ++gp)
        A[g][gp] = -S[g][gp];
      for (unsigned int gp = g + 1; gp < num_groups_; ++gp)
        B[g][gp] = S[g][gp];
    }
    MatDbl Ainv = chi_math::Inverse(A);
    // finally, obtain the iteration matrix
    MatDbl C_ = chi_math::MatMul(Ainv, B);
    // Perform power iteration
    VecDbl E(num_groups_, 1.0);
    double rho = chi_math::PowerIteration(C_, E, 10000, 1.0e-12);

    // Compute two-grid diffusion quantities
    // normalize spectrum
    std::vector<double> spectrum(num_groups_, 1.0);
    double sum = 0.0;
    for (unsigned int g = 0; g < num_groups_; ++g)
      sum += std::fabs(E[g]);
    for (unsigned int g = 0; g < num_groups_; ++g)
      spectrum[g] = std::fabs(E[g]) / sum;
    // D ave and Sigma_a ave
    double collapsed_D = 0.0;
    double collapsed_sig_a = 0.0;
    for (unsigned int g = last_fast_group_; g < num_groups_; ++g)
    {
      collapsed_D += diffusion_coeff[g] * spectrum[g];
      collapsed_sig_a += sigma_t[g] * spectrum[g];
      for (unsigned int gp = last_fast_group_; gp < num_groups_; ++gp)
        collapsed_sig_a -= S[g][gp] * spectrum[gp];
    }
    // Verbose output the spectrum
    Chi::log.Log0Verbose1() << "Fundamental eigen-value: " << rho;
    std::stringstream outstr;
    for (auto &xi: spectrum)
      outstr << xi << '\n';
    Chi::log.Log0Verbose1() << outstr.str();  // jcr verbose1

//    std::stringstream outstr2;
//    for (auto &xi: diffusion_coeff)
//      outstr2 << xi << '\n';
//    chi::log.Log0Verbose0() << outstr2.str();  // jcr verbose1
//
//    std::cout << "collapsed = " << collapsed_sig_a
//    <<", "<< collapsed_D << std::endl;
//    chi::Exit(12345);

    const auto mat_id = mat_id_xs.first;
    map_mat_id_2_tginfo.insert(
      std::make_pair(mat_id, TwoGridCollapsedInfo{collapsed_D,
                                                  collapsed_sig_a,
                                                  spectrum} ) );

  }// end loop over materials

}
