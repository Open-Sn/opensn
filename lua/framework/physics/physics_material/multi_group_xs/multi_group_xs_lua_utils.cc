#include "framework/lua.h"
#include <iostream>

#include "framework/physics/physics_namespace.h"
#include "framework/physics/physics_material/multi_group_xs/single_state_mgxs.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "multi_group_xs_lua_utils.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiPhysicsTransportXSCreate);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSSet);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSMakeCombined);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSSetCombined);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSGet);
RegisterLuaFunctionAsIs(chiPhysicsTransportXSExportToChiTechFormat);

RegisterLuaConstantAsIs(SINGLE_VALUE, Varying(0));
RegisterLuaConstantAsIs(FROM_ARRAY, Varying(1));
RegisterLuaConstantAsIs(SIMPLEXS0, Varying(20));
RegisterLuaConstantAsIs(SIMPLEXS1, Varying(21));
RegisterLuaConstantAsIs(EXISTING, Varying(22));
RegisterLuaConstantAsIs(CHI_XSFILE, Varying(23));

namespace
{

void
MultiGroupXSPushLuaTable(lua_State* L, std::shared_ptr<MultiGroupXS> xs)
{
  // General data
  lua_newtable(L);
  lua_pushstring(L, "is_empty");
  lua_pushboolean(L, false);
  lua_settable(L, -3);

  lua_pushstring(L, "num_groups");
  lua_pushinteger(L, static_cast<lua_Integer>(xs->NumGroups()));
  lua_settable(L, -3);

  lua_pushstring(L, "scattering_order");
  lua_pushinteger(L, static_cast<lua_Integer>(xs->ScatteringOrder()));
  lua_settable(L, -3);

  lua_pushstring(L, "num_precursors");
  lua_pushinteger(L, static_cast<lua_Integer>(xs->NumPrecursors()));
  lua_settable(L, -3);

  lua_pushstring(L, "is_fissionable");
  lua_pushboolean(L, xs->IsFissionable());
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

  Push1DXS(xs->SigmaTotal(), "sigma_t");
  Push1DXS(xs->SigmaAbsorption(), "sigma_a");
  Push1DXS(xs->SigmaFission(), "sigma_f");
  Push1DXS(xs->NuSigmaF(), "nu_sigma_f");
  Push1DXS(xs->NuPromptSigmaF(), "nu_prompt_sigma_f");
  Push1DXS(xs->NuDelayedSigmaF(), "nu_delayed_sigma_f");
  Push1DXS(xs->InverseVelocity(), "inv_velocity");

  // Emission spectra
  std::vector<std::vector<double>> chi_delayed;
  for (unsigned int g = 0; g < xs->NumGroups(); ++g)
  {
    std::vector<double> vals;
    for (const auto& precursor : xs->Precursors())
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
    for (const auto& precursor : xs->Precursors())
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
    for (const auto& precursor : xs->Precursors())
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
    for (const auto& matrix : xs->TransferMatrices())
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
    for (const auto& prod : xs->ProductionMatrix())
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
  Push1DXS(xs->DiffusionCoefficient(), "diffusion_coeff");
  Push1DXS(xs->SigmaRemoval(), "sigma_removal");
  Push1DXS(xs->SigmaSGtoG(), "sigma_s_gtog");
}

} // namespace

int
chiPhysicsTransportXSCreate(lua_State* L)
{
  auto xs = std::make_shared<SingleStateMGXS>();
  opensn::Chi::multigroup_xs_stack.push_back(xs);

  const size_t index = opensn::Chi::multigroup_xs_stack.size() - 1;
  lua_pushinteger(L, static_cast<lua_Integer>(index));
  return 1;
}

int
chiPhysicsTransportXSSet(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  const int num_args = lua_gettop(L);

  if (num_args < 3) LuaPostArgAmountError(fname, 3, num_args);

  // Process xs handle
  LuaCheckIntegerValue(fname, L, 1);
  const int handle = lua_tointeger(L, 1);

  // Process operation id
  LuaCheckIntegerValue(fname, L, 2);
  const int operation_index = lua_tointeger(L, 2);

  std::shared_ptr<SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<SingleStateMGXS>(
      opensn::Chi::GetStackItemPtr(opensn::Chi::multigroup_xs_stack, handle));
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }

  // Process operation
  using OpType = OperationType;
  if (operation_index == static_cast<int>(OpType::SIMPLEXS0))
  {
    if (num_args != 4) LuaPostArgAmountError(fname, 4, num_args);

    const int n_grps = lua_tointeger(L, 3);
    const double sigma_t = lua_tonumber(L, 4);

    xs->MakeSimple0(n_grps, sigma_t);
  }
  else if (operation_index == static_cast<int>(OpType::SIMPLEXS1))
  {
    if (num_args != 5) LuaPostArgAmountError(fname, 5, num_args);

    const int n_grps = lua_tointeger(L, 3);
    const double sigma_t = lua_tonumber(L, 4);
    const double c = lua_tonumber(L, 5);

    xs->MakeSimple1(n_grps, sigma_t, c);
  }
  else if (operation_index == static_cast<int>(OpType::CHI_XSFILE))
  {
    if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

    const char* file_name_c = lua_tostring(L, 3);

    xs->MakeFromChiXSFile(std::string(file_name_c));
  }
  else
  {
    opensn::log.LogAllError() << "Unsupported operation in " << fname << ". " << operation_index;
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
}

int
chiPhysicsTransportXSGet(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  const int num_args = lua_gettop(L);

  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  // Process xs handle
  LuaCheckNilValue(fname, L, 1);
  const int handle = lua_tointeger(L, 1);

  std::shared_ptr<SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<SingleStateMGXS>(
      opensn::Chi::GetStackItemPtr(opensn::Chi::multigroup_xs_stack, handle));
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }

  MultiGroupXSPushLuaTable(L, xs);

  return 1;
}

int
chiPhysicsTransportXSMakeCombined(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  // Process table
  LuaCheckTableValue(fname, L, 1);
  std::vector<std::pair<int, double>> combinations;
  for (int v = 0; v < lua_rawlen(L, 1); ++v)
  {
    lua_pushnumber(L, v + 1);
    lua_gettable(L, 1);

    ChiInvalidArgumentIf(not lua_istable(L, -1),
                         "The elements of the supplied Lua array must themselves also be "
                         "Lua array containing a xs handle and scalar multiplier.");

    // Process xs handle
    lua_pushinteger(L, 1);
    lua_gettable(L, -2);
    LuaCheckIntegerValue(fname + "A1:E1", L, -1);
    const int handle = lua_tointeger(L, -1);
    lua_pop(L, 1);

    // Process scalar multiplier
    lua_pushinteger(L, 2);
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
  auto new_xs = std::make_shared<SingleStateMGXS>();

  new_xs->MakeCombined(combinations);

  opensn::Chi::multigroup_xs_stack.push_back(new_xs);
  const lua_Integer num_xs = opensn::Chi::multigroup_xs_stack.size();
  lua_pushinteger(L, num_xs - 1);

  return 1;
}

int
chiPhysicsTransportXSSetCombined(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  int num_args = lua_gettop(L);

  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  // Process xs handle
  LuaCheckIntegerValue(fname, L, 1);
  const int xs_handle = lua_tointeger(L, 1);

  std::shared_ptr<SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<SingleStateMGXS>(
      opensn::Chi::GetStackItemPtr(opensn::Chi::multigroup_xs_stack, xs_handle));
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }

  // Process table
  LuaCheckTableValue(fname, L, 2);
  std::vector<std::pair<int, double>> combinations;
  for (int v = 0; v < lua_rawlen(L, 2); ++v)
  {
    lua_pushnumber(L, v + 1);
    lua_gettable(L, 1);

    ChiInvalidArgumentIf(not lua_istable(L, -1),
                         "The elements of the supplied Lua array must themselves also be "
                         "Lua array containing a xs handle and scalar multiplier.");

    // Process xs handle
    lua_pushinteger(L, 1);
    lua_gettable(L, -2);
    LuaCheckIntegerValue(fname + ":A1:E1", L, -1);
    const int handle = lua_tonumber(L, -1);
    lua_pop(L, 1);

    // Process scalar multiplier
    lua_pushinteger(L, 2);
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

  xs->MakeCombined(combinations);

  return 0;
}

int
chiPhysicsTransportXSExportToChiTechFormat(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  const int num_args = lua_gettop(L);

  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  // Process xs handle
  LuaCheckIntegerValue(fname, L, 1);
  const int handle = lua_tointeger(L, 1);

  std::shared_ptr<MultiGroupXS> xs;
  try
  {
    xs = opensn::Chi::GetStackItemPtr(opensn::Chi::multigroup_xs_stack, handle);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }

  // Process file name
  LuaCheckNilValue(fname, L, 2);
  std::string file_name = lua_tostring(L, 2);

  xs->ExportToChiXSFile(file_name);

  return 0;
}
