// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework//materials/multi_group_xs/multi_group_xs_lua_utils.h"
#include "lua/framework/console/console.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/material.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iostream>

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(XSCreate, xs, Create);
RegisterLuaFunctionInNamespace(XSSet, xs, Set);
RegisterLuaFunctionInNamespace(XSMakeCombined, xs, MakeCombined);
RegisterLuaFunctionInNamespace(XSSetCombined, xs, SetCombined);
RegisterLuaFunctionInNamespace(XSMakeScaled, xs, MakeScaled);
RegisterLuaFunctionInNamespace(XSSetScalingFactor, xs, SetScalingFactor);
RegisterLuaFunctionInNamespace(XSGet, xs, Get);
RegisterLuaFunctionInNamespace(XSExportToOpenSnFormat, xs, ExportToOpenSnFormat);

RegisterLuaConstant(SINGLE_VALUE, Varying(0));
RegisterLuaConstant(FROM_ARRAY, Varying(1));
RegisterLuaConstant(SIMPLE_ONE_GROUP, Varying(20));
RegisterLuaConstant(EXISTING, Varying(22));
RegisterLuaConstant(OPENSN_XSFILE, Varying(23));
RegisterLuaConstant(OPENMC_XSLIB, Varying(24));

namespace
{

void
MultiGroupXSPushLuaTable(lua_State* L, std::shared_ptr<MultiGroupXS> xs)
{
  // General data
  lua_newtable(L);
  LuaPushTableKey(L, "is_empty", false);
  LuaPushTableKey(L, "num_groups", xs->NumGroups());
  LuaPushTableKey(L, "scattering_order", xs->ScatteringOrder());
  LuaPushTableKey(L, "num_precursors", xs->NumPrecursors());
  LuaPushTableKey(L, "is_fissionable", xs->IsFissionable());
  LuaPushTableKey(L, "scaling_factor", xs->ScalingFactor());
  LuaPushTableKey(L, "sigma_t", xs->SigmaTotal());
  LuaPushTableKey(L, "sigma_a", xs->SigmaAbsorption());
  LuaPushTableKey(L, "sigma_f", xs->SigmaFission());
  LuaPushTableKey(L, "chi", xs->Chi());
  LuaPushTableKey(L, "nu_sigma_f", xs->NuSigmaF());
  LuaPushTableKey(L, "nu_prompt_sigma_f", xs->NuPromptSigmaF());
  LuaPushTableKey(L, "nu_delayed_sigma_f", xs->NuDelayedSigmaF());
  LuaPushTableKey(L, "inv_velocity", xs->InverseVelocity());

  // Emission spectra
  std::vector<std::vector<double>> chi_delayed;
  for (unsigned int g = 0; g < xs->NumGroups(); ++g)
  {
    std::vector<double> vals;
    for (const auto& precursor : xs->Precursors())
      vals.push_back(precursor.emission_spectrum[g]);
    chi_delayed.push_back(vals);
  }

  LuaPushTableKey(L, "chi_delayed", chi_delayed);

  // Precursor data
  LuaPush(L, "precursor_decay_constants");
  lua_newtable(L);
  {
    unsigned int j = 0;
    for (const auto& precursor : xs->Precursors())
      LuaPushTableKey(L, ++j, precursor.decay_constant);
  }
  lua_settable(L, -3);

  LuaPush(L, "precursor_fractional_yields");
  lua_newtable(L);
  {
    unsigned int j = 0;
    for (const auto& precursor : xs->Precursors())
      LuaPushTableKey(L, ++j, precursor.fractional_yield);
  }
  lua_settable(L, -3);

  // Transfer matrices
  LuaPush(L, "transfer_matrix");
  lua_newtable(L);
  {
    unsigned int ell = 0;
    for (const auto& matrix : xs->TransferMatrices())
    {
      LuaPush(L, ++ell);
      lua_newtable(L);
      {
        for (unsigned int g = 0; g < matrix.NumRows(); ++g)
        {
          const auto& col_indices = matrix.rowI_indices_[g];
          const auto& col_values = matrix.rowI_values_[g];
          size_t num_vals = col_values.size();

          LuaPush(L, g + 1);
          lua_newtable(L);
          {
            for (unsigned int gg = 0; gg < num_vals; ++gg)
              LuaPushTableKey(L, col_indices[gg] + 1, col_values[gg]);
          }
          lua_settable(L, -3);
        }
      }
      lua_settable(L, -3);
    }
  }
  lua_settable(L, -3);

  // Production matrix
  LuaPushTableKey(L, "production_matrix", xs->ProductionMatrix());

  // Push diffusion quantities
  LuaPushTableKey(L, "diffusion_coeff", xs->DiffusionCoefficient());
  LuaPushTableKey(L, "sigma_removal", xs->SigmaRemoval());
  LuaPushTableKey(L, "sigma_s_gtog", xs->SigmaSGtoG());
}

} // namespace

int
XSCreate(lua_State* L)
{
  auto xs = std::make_shared<MultiGroupXS>();
  opensn::multigroup_xs_stack.push_back(xs);

  const size_t index = opensn::multigroup_xs_stack.size() - 1;
  return LuaReturn(L, index);
}

int
XSSet(lua_State* L)
{
  const std::string fname = "xs.Set";
  LuaCheckArgs<int, int>(L, fname);

  // Process xs handle
  const auto handle = LuaArg<int>(L, 1);
  // Process operation id
  const auto operation_index = LuaArg<int>(L, 2);

  std::shared_ptr<MultiGroupXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<MultiGroupXS>(
      opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle));
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross-section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }

  // Process operation
  if (operation_index == static_cast<int>(OperationType::SIMPLE_ONE_GROUP))
  {
    LuaCheckArgs<int, int, double, double>(L, fname);

    const auto sigma_t = LuaArg<double>(L, 3);
    const auto c = LuaArg<double>(L, 4);

    xs->Initialize(sigma_t, c);
  }
  else if (operation_index == static_cast<int>(OperationType::OPENSN_XSFILE))
  {
    LuaCheckArgs<int, int, std::string>(L, fname);

    auto file_name = LuaArg<std::string>(L, 3);

    xs->Initialize(file_name);
  }
  else if (operation_index == static_cast<int>(OperationType::OPENMC_XSLIB))
  {
    LuaCheckArgs<int, int, std::string, std::string>(L, fname);
    auto file_name = LuaArg<std::string>(L, 3);
    auto temperature = LuaArg<double>(L, 4);
    const auto xs_data_name = LuaArgOptional<std::string>(L, 5, "set1");
    xs->Initialize(file_name, xs_data_name, temperature);
  }
  else
  {
    opensn::log.LogAllError() << "Unsupported operation in " << fname << ". " << operation_index;
    opensn::Exit(EXIT_FAILURE);
  }
  return LuaReturn(L);
}

int
XSGet(lua_State* L)
{
  const std::string fname = "xs.Get";
  LuaCheckArgs<int>(L, fname);

  // Process xs handle
  const auto handle = LuaArg<int>(L, 1);

  std::shared_ptr<MultiGroupXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<MultiGroupXS>(
      opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle));
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross-section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }

  MultiGroupXSPushLuaTable(L, xs);

  return 1;
}

int
XSMakeCombined(lua_State* L)
{
  const std::string fname = "xs.MakeCombined";

  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, L, 1, num_args);

  // Process table
  LuaCheckTableValue(fname, L, 1);
  std::vector<std::pair<int, double>> combinations;
  for (int v = 0; v < lua_rawlen(L, 1); ++v)
  {
    LuaPush(L, v + 1);
    lua_gettable(L, 1);

    OpenSnInvalidArgumentIf(not lua_istable(L, -1),
                            "The elements of the supplied Lua array must themselves also be "
                            "Lua array containing a xs handle and scalar multiplier.");

    // Process xs handle
    LuaPush(L, 1);
    lua_gettable(L, -2);
    LuaCheckIntegerValue(fname + "A1:E1", L, -1);
    const int handle = lua_tointeger(L, -1);
    lua_pop(L, 1);

    // Process scalar multiplier
    LuaPush(L, 2);
    lua_gettable(L, -2);
    LuaCheckNumberValue(fname + ":A1:E2", L, -1);
    const double scalar = lua_tonumber(L, -1);
    lua_pop(L, 1);

    combinations.emplace_back(handle, scalar);
    lua_pop(L, 1); // pop off table
  }

  // Print out table
  opensn::log.Log() << "Generating XS with following combination:";
  for (auto& elem : combinations)
    opensn::log.Log() << " Element handle: " << elem.first << " scalar value: " << elem.second;

  // Make the new cross section
  auto new_xs = std::make_shared<MultiGroupXS>();

  new_xs->Initialize(combinations);

  opensn::multigroup_xs_stack.push_back(new_xs);
  auto num_xs = opensn::multigroup_xs_stack.size();
  return LuaReturn(L, num_xs - 1);
}

int
XSSetCombined(lua_State* L)
{
  const std::string fname = "xs.SetCombined";

  int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(fname, L, 2, num_args);

  // Process xs handle
  const auto xs_handle = LuaArg<int>(L, 1);

  std::shared_ptr<MultiGroupXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<MultiGroupXS>(
      opensn::GetStackItemPtr(opensn::multigroup_xs_stack, xs_handle));
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross-section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }

  // Process table
  LuaCheckTableValue(fname, L, 2);
  std::vector<std::pair<int, double>> combinations;
  for (int v = 0; v < lua_rawlen(L, 2); ++v)
  {
    LuaPush(L, v + 1);
    lua_gettable(L, 1);

    OpenSnInvalidArgumentIf(not lua_istable(L, -1),
                            "The elements of the supplied Lua array must themselves also be "
                            "Lua array containing a xs handle and scalar multiplier.");

    // Process xs handle
    LuaPush(L, 1);
    lua_gettable(L, -2);
    LuaCheckIntegerValue(fname + ":A1:E1", L, -1);
    const int handle = lua_tonumber(L, -1);
    lua_pop(L, 1);

    // Process scalar multiplier
    LuaPush(L, 2);
    lua_gettable(L, -2);
    LuaCheckNumberValue(fname + ":A1:E2", L, -1);
    const double scalar = lua_tonumber(L, -1);
    lua_pop(L, 1);

    combinations.emplace_back(handle, scalar);
    lua_pop(L, 1); // pop off table
  }

  // Print out table
  opensn::log.Log() << "Setting XS with following combination:";
  for (auto& elem : combinations)
    opensn::log.Log() << "  Element handle: " << elem.first << " scalar value: " << elem.second;

  xs->Initialize(combinations);

  return LuaReturn(L);
}

int
XSMakeScaled(lua_State* L)
{
  const std::string fname = "xs.MakeScaled";
  LuaCheckArgs<int, double>(L, fname);

  const auto handle = LuaArg<int>(L, 1);
  auto factor = LuaArg<double>(L, 2);

  std::shared_ptr<MultiGroupXS> xs;
  try
  {
    xs = opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle);
  }
  catch (const std::out_of_range& o)
  {
    OpenSnInvalidArgument("Invalid handle for cross sections in call to " + fname + ".");
  }

  auto new_xs = std::shared_ptr<MultiGroupXS>(xs);
  new_xs->SetScalingFactor(factor);
  opensn::multigroup_xs_stack.push_back(new_xs);
  auto num_xs = opensn::multigroup_xs_stack.size();
  return LuaReturn(L, num_xs - 1);
}

int
XSSetScalingFactor(lua_State* L)
{
  const std::string fname = "xs.SetScalingFactor";
  LuaCheckArgs<int, double>(L, fname);

  const auto handle = LuaArg<int>(L, 1);
  auto factor = LuaArg<double>(L, 2);

  std::shared_ptr<MultiGroupXS> xs;
  try
  {
    xs = opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle);
  }
  catch (const std::out_of_range& o)
  {
    OpenSnInvalidArgument("Invalid handle for cross sections in call to " + fname + ".");
  }

  xs->SetScalingFactor(factor);
  return LuaReturn(L);
}

int
XSExportToOpenSnFormat(lua_State* L)
{
  const std::string fname = "xs.ExportToOpenSnFormat";
  LuaCheckArgs<int, std::string>(L, fname);

  const auto handle = LuaArg<int>(L, 1);
  auto file_name = LuaArg<std::string>(L, 2);

  std::shared_ptr<MultiGroupXS> xs;
  try
  {
    xs = opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross-section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }
  xs->ExportToOpenSnXSFile(file_name);

  return LuaReturn(L);
}

} // namespace opensnlua
