#include "multigroup_xs.h"

#include "chi_log.h"

//###################################################################
/**Pushes all of the relevant items of the transport xs to a lua table.*/
void
chi_physics::MultiGroupXS::PushLuaTable(lua_State* L) const
{
  //================================================== General data
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

  //================================================== 1D cross sections
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

  //================================================== Emission spectra
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

  //================================================== Precursor data
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

  //================================================== Transfer matrices
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

  //================================================== Production matrix
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

  //================================================== Push diffusion quantities
  Push1DXS(DiffusionCoefficient(), "diffusion_coeff");
  Push1DXS(SigmaRemoval(), "sigma_removal");
  Push1DXS(SigmaSGtoG(), "sigma_s_gtog");
}
