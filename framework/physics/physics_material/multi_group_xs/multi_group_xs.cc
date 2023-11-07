#include "framework/physics/physics_material/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include <iostream>

namespace chi_physics
{

void
MultiGroupXS::ExportToChiXSFile(const std::string& file_name, const double fission_scaling) const
{
  Chi::log.Log() << "Exporting transport cross section to file: " << file_name;

  // Define utility functions

  /**Lambda to print a 1D-xs*/
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

      if (not proceed) return;
    }

    ofile << "\n";
    ofile << prefix << "_BEGIN\n";
    {
      unsigned int g = 0;
      for (auto val : xs)
        ofile << g++ << " " << val << "\n";
    }
    ofile << prefix << "_END\n";
  };

  // Open the output file
  std::ofstream ofile(file_name);

  // Write the header info
  std::vector<double> nu, nu_prompt, nu_delayed;
  const auto& sigma_f = SigmaFission();
  const auto& nu_sigma_f = NuSigmaF();
  const auto& nu_prompt_sigma_f = NuPromptSigmaF();
  const auto& nu_delayed_sigma_f = NuDelayedSigmaF();
  for (unsigned int g = 0; g < NumGroups(); ++g)
  {
    if (NumPrecursors() > 0)
    {
      nu_prompt.push_back(nu_prompt_sigma_f[g] / sigma_f[g]);
      nu_delayed.push_back(nu_delayed_sigma_f[g] / sigma_f[g]);
    }
    nu.push_back(nu_sigma_f[g] / sigma_f[g]);
  }

  std::vector<double> decay_constants, fractional_yields;
  for (const auto& precursor : Precursors())
  {
    decay_constants.push_back(precursor.decay_constant);
    fractional_yields.push_back(precursor.fractional_yield);
  }

  ofile << "# Exported cross section from ChiTech\n";
  ofile << "# Date: " << chi::Timer::GetLocalDateTimeString() << "\n";
  ofile << "NUM_GROUPS " << NumGroups() << "\n";
  ofile << "NUM_MOMENTS " << ScatteringOrder() + 1 << "\n";
  if (NumPrecursors() > 0) ofile << "NUM_PRECURSORS " << NumPrecursors() << "\n";

  // basic cross section data
  Print1DXS(ofile, "SIGMA_T", SigmaTotal(), 1.0e-20);
  Print1DXS(ofile, "SIGMA_A", SigmaAbsorption(), 1.0e-20);

  // fission data
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
      Print1DXS(ofile, "NU_PROMPT", nu_prompt, 1.0e-20);
      Print1DXS(ofile, "NU_DELAYED", nu_delayed, 1.0e-20);
      //      Print1DXS(ofile, "CHI_PROMPT", chi_prompt, 1.0e-20);

      ofile << "\nCHI_DELAYED_BEGIN\n";
      const auto& precursors = Precursors();
      for (unsigned int j = 0; j < NumPrecursors(); ++j)
        for (unsigned int g = 0; g < NumGroups(); ++g)
          ofile << "G_PRECURSOR_VAL"
                << " " << g << " " << j << " " << precursors[j].emission_spectrum[g] << "\n";
      ofile << "CHI_DELAYED_END\n";

      Print1DXS(ofile, "PRECURSOR_DECAY_CONSTANTS", decay_constants, 1.0e-20);
      Print1DXS(ofile, "PRECURSOR_FRACTIONAL_YIELDS", fractional_yields, 1.0e-20);
    }
    else
    {
      Print1DXS(ofile, "NU", nu, 1.0e-20);
      //      Print1DXS(ofile, "CHI", chi, 1.0e-20);
    }
  }

  // inverse speed data
  if (not InverseVelocity().empty()) Print1DXS(ofile, "INV_VELOCITY", InverseVelocity(), 1.0e-20);

  // transfer matrices
  if (not TransferMatrices().empty())
  {
    ofile << "\n";
    ofile << "TRANSFER_MOMENTS_BEGIN\n";
    for (size_t ell = 0; ell < TransferMatrices().size(); ++ell)
    {
      if (ell == 0) ofile << "#Zeroth moment (l=0)\n";
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

  if (not ProductionMatrix().empty())
  {
    const auto& F = ProductionMatrix();

    ofile << "\n";
    ofile << "PRODUCTION_MATRIX_BEGIN\n";
    for (unsigned int g = 0; g < NumGroups(); ++g)
    {
      const auto& prod = F[g];
      for (unsigned int gp = 0; gp < NumGroups(); ++gp)
      {
        const double value = fission_scaling != 1.0 ? prod[gp] * fission_scaling : prod[gp];

        ofile << "G_GPRIME_VAL " << g << " " << gp << " " << value << "\n";
      }
    }
  }

  ofile.close();

  Chi::log.Log0Verbose1() << "Done exporting transport "
                             "cross section to file: "
                          << file_name;
}

#ifdef OPENSN_WITH_LUA
void
MultiGroupXS::PushLuaTable(lua_State* L) const
{
  // General data
  lua_newtable(L);
  lua_pushstring(L, "is_empty");
  lua_pushboolean(L, false);
  lua_settable(L, -3);

  lua_pushstring(L, "num_groups");
  lua_pushinteger(L, static_cast<lua_Integer>(NumGroups()));
  lua_settable(L, -3);

  lua_pushstring(L, "scattering_order");
  lua_pushinteger(L, static_cast<lua_Integer>(ScatteringOrder()));
  lua_settable(L, -3);

  lua_pushstring(L, "num_precursors");
  lua_pushinteger(L, static_cast<lua_Integer>(NumPrecursors()));
  lua_settable(L, -3);

  lua_pushstring(L, "is_fissionable");
  lua_pushboolean(L, IsFissionable());
  lua_settable(L, -3);

  // 1D cross sections
  auto Push1DXS = [L](const std::vector<double>& xs, const std::string& name)
  {
    lua_pushstring(L, name.c_str());
    lua_newtable(L);
    {
      unsigned int g = 0;
      for (const auto& val : xs)
      {
        lua_pushinteger(L, ++g);
        lua_pushnumber(L, val);
        lua_settable(L, -3);
      }
    }
    lua_settable(L, -3);
  };

  Push1DXS(SigmaTotal(), "sigma_t");
  Push1DXS(SigmaAbsorption(), "sigma_a");
  Push1DXS(SigmaFission(), "sigma_f");
  Push1DXS(NuSigmaF(), "nu_sigma_f");
  Push1DXS(NuPromptSigmaF(), "nu_prompt_sigma_f");
  Push1DXS(NuDelayedSigmaF(), "nu_delayed_sigma_f");
  Push1DXS(InverseVelocity(), "inv_velocity");

  // Emission spectra
  std::vector<std::vector<double>> chi_delayed;
  for (unsigned int g = 0; g < NumGroups(); ++g)
  {
    std::vector<double> vals;
    for (const auto& precursor : Precursors())
      vals.push_back(precursor.emission_spectrum[g]);
    chi_delayed.push_back(vals);
  }

  lua_pushstring(L, "chi_delayed");
  lua_newtable(L);
  {
    unsigned int g = 0;
    for (const auto& emission_g : chi_delayed)
    {
      lua_pushinteger(L, ++g);
      lua_newtable(L);
      {
        unsigned int j = 0;
        for (const auto& val : emission_g)
        {
          lua_pushinteger(L, ++j);
          lua_pushnumber(L, val);
          lua_settable(L, -3);
        }
      }
      lua_settable(L, -3);
    }
  }
  lua_settable(L, -3);

  // Precursor data
  lua_pushstring(L, "precursor_decay_constants");
  lua_newtable(L);
  {
    unsigned int j = 0;
    for (const auto& precursor : Precursors())
    {
      lua_pushinteger(L, ++j);
      lua_pushnumber(L, precursor.decay_constant);
      lua_settable(L, -3);
    }
  }
  lua_settable(L, -3);

  lua_pushstring(L, "precursor_fractional_yields");
  lua_newtable(L);
  {
    unsigned int j = 0;
    for (const auto& precursor : Precursors())
    {
      lua_pushinteger(L, ++j);
      lua_pushnumber(L, precursor.fractional_yield);
      lua_settable(L, -3);
    }
  }
  lua_settable(L, -3);

  // Transfer matrices
  lua_pushstring(L, "transfer_matrix");
  lua_newtable(L);
  {
    unsigned int ell = 0;
    for (const auto& matrix : TransferMatrices())
    {
      lua_pushinteger(L, ++ell);
      lua_newtable(L);
      {
        for (unsigned int g = 0; g < matrix.NumRows(); ++g)
        {
          const auto& col_indices = matrix.rowI_indices_[g];
          const auto& col_values = matrix.rowI_values_[g];
          size_t num_vals = col_values.size();

          lua_pushinteger(L, g + 1);
          lua_newtable(L);
          {
            for (unsigned int gg = 0; gg < num_vals; ++gg)
            {
              lua_pushinteger(L, static_cast<long long>(col_indices[gg]) + 1);
              lua_pushnumber(L, col_values[gg]);
              lua_settable(L, -3);
            }
            lua_settable(L, -3);
          }
        }
      }
      lua_settable(L, -3);
    }
  }
  lua_settable(L, -3);

  // Production matrix
  lua_pushstring(L, "production_matrix");
  lua_newtable(L);
  {
    unsigned int g = 0;
    for (const auto& prod : ProductionMatrix())
    {
      lua_pushinteger(L, ++g);
      lua_newtable(L);
      {
        unsigned int gp = 0;
        for (const auto& val : prod)
        {
          lua_pushinteger(L, ++gp);
          lua_pushnumber(L, val);
          lua_settable(L, -3);
        }
        lua_settable(L, -3);
      }
    } // for g
  }
  lua_settable(L, -3);

  // Push diffusion quantities
  Push1DXS(DiffusionCoefficient(), "diffusion_coeff");
  Push1DXS(SigmaRemoval(), "sigma_removal");
  Push1DXS(SigmaSGtoG(), "sigma_s_gtog");
}
#endif

} // namespace chi_physics
