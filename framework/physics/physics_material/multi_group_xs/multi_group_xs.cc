#include "framework/physics/physics_material/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include <iostream>

namespace opensn
{

void
MultiGroupXS::ExportToOpenSnXSFile(const std::string& file_name,
                                   const double fission_scaling /* = 1.0 */) const
{
  log.Log() << "Exporting transport cross section to file: " << file_name;

  // Lambda to print a 1D cross section
  auto Print1DXS = [](std::ofstream& ofile,
                      const std::string& prefix,
                      const std::vector<double>& xs,
                      double min_value = -1.0)
  {
    bool proceed = false;
    if (min_value >= 0.0)
    {
      for (auto val : xs)
        if (val > min_value)
        {
          proceed = true;
          break;
        }

      if (not proceed)
        return;
    }

    ofile << "\n";
    ofile << prefix << "_BEGIN\n";
    {
      size_t g = 0;
      for (auto val : xs)
        ofile << g++ << " " << val << "\n";
    }
    ofile << prefix << "_END\n";
  };

  // Open the output file
  std::ofstream ofile(file_name);

  // Write the header info
  ofile << "# Exported cross section from OpenSn\n";
  ofile << "# Date: " << Timer::GetLocalDateTimeString() << "\n";
  ofile << "NUM_GROUPS " << NumGroups() << "\n";
  ofile << "NUM_MOMENTS " << ScatteringOrder() + 1 << "\n";
  if (NumPrecursors() > 0)
    ofile << "NUM_PRECURSORS " << NumPrecursors() << "\n";

  // Basic cross section data
  Print1DXS(ofile, "SIGMA_T", SigmaTotal(), 1.0e-20);
  Print1DXS(ofile, "SIGMA_A", SigmaAbsorption(), 1.0e-20);

  // Fission data
  if (not SigmaFission().empty())
  {
    std::vector<double> scaled_sigma_f = SigmaFission();
    if (fission_scaling != 1.0)
    {
      for (auto& val : scaled_sigma_f)
        val *= fission_scaling;
    }

    Print1DXS(ofile, "SIGMA_F", scaled_sigma_f, 1.0e-20);
    if (NumPrecursors() > 0)
    {
      // Compute prompt/delayed fission neutron yields
      const auto& sigma_f = SigmaFission();
      const auto& nu_prompt_sigma_f = NuPromptSigmaF();
      const auto& nu_delayed_sigma_f = NuDelayedSigmaF();

      std::vector<double> nu_prompt, nu_delayed;
      for (size_t g = 0; g < NumGroups(); ++g)
      {
        nu_prompt.emplace_back(sigma_f[g] > 0.0 ? nu_prompt_sigma_f[g] / sigma_f[g] : 0.0);
        nu_delayed.emplace_back(sigma_f[g] > 0.0 ? nu_delayed_sigma_f[g] / sigma_f[g] : 0.0);
      }

      // Get decay constants and fractional yields
      std::vector<double> lambda, gamma;
      for (const auto& precursor : Precursors())
      {
        lambda.emplace_back(precursor.decay_constant);
        gamma.emplace_back(precursor.fractional_yield);
      }

      Print1DXS(ofile, "NU_PROMPT", nu_prompt, 1.0e-20);
      Print1DXS(ofile, "NU_DELAYED", nu_delayed, 1.0e-20);

      ofile << "\nCHI_DELAYED_BEGIN\n";
      const auto& precursors = Precursors();
      for (size_t j = 0; j < NumPrecursors(); ++j)
        for (size_t g = 0; g < NumGroups(); ++g)
          ofile << "G_PRECURSOR_VAL"
                << " " << g << " " << j << " " << precursors[j].emission_spectrum[g] << "\n";
      ofile << "CHI_DELAYED_END\n";

      Print1DXS(ofile, "PRECURSOR_DECAY_CONSTANTS", lambda, 1.0e-20);
      Print1DXS(ofile, "PRECURSOR_FRACTIONAL_YIELDS", gamma, 1.0e-20);
    }
    else
    {
      // Compute the average total fission neutron yield
      const auto& sigma_f = SigmaFission();
      const auto& nu_sigma_f = NuSigmaF();

      std::vector<double> nu;
      for (size_t g = 0; g < NumGroups(); ++g)
        nu.emplace_back(sigma_f[g] > 0.0 ? nu_sigma_f[g] / sigma_f[g] : 0.0);

      Print1DXS(ofile, "NU", nu, 1.0e-20);
    }
  }

  // Inverse velocity data
  if (not InverseVelocity().empty())
    Print1DXS(ofile, "INV_VELOCITY", InverseVelocity(), 1.0e-20);

  // Transfer matrices
  if (not TransferMatrices().empty())
  {
    ofile << "\n";
    ofile << "TRANSFER_MOMENTS_BEGIN\n";
    for (size_t ell = 0; ell < TransferMatrices().size(); ++ell)
    {
      if (ell == 0)
        ofile << "#Zeroth moment (l=0)\n";
      else
        ofile << "#(l=" << ell << ")\n";

      const auto& matrix = TransferMatrix(ell);

      for (size_t g = 0; g < matrix.rowI_values_.size(); ++g)
      {
        const auto& col_indices = matrix.rowI_indices_[g];
        const auto& col_values = matrix.rowI_values_[g];
        for (size_t k = 0; k < col_indices.size(); ++k)
          ofile << "M_GPRIME_G_VAL " << ell << " " << col_indices[k] << " " << g << " "
                << col_values[k] << "\n";
      } // for g

      ofile << "\n";
    } // for ell
    ofile << "TRANSFER_MOMENTS_END\n";
  } // if has transfer matrices

  // Write production matrix
  if (not ProductionMatrix().empty())
  {
    ofile << "\n";
    ofile << "PRODUCTION_MATRIX_BEGIN\n";
    for (size_t g = 0; g < NumGroups(); ++g)
      for (size_t gp = 0; gp < NumGroups(); ++gp)
        ofile << "G_GPRIME_VAL " << g << " " << gp << " "
              << fission_scaling * ProductionMatrix()[g][gp] << "\n";
  }
  ofile.close();

  log.Log0Verbose1() << "Done exporting transport cross section to file: " << file_name;
}

} // namespace opensn
