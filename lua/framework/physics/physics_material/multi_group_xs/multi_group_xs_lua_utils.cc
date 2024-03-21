#include "multi_group_xs_lua_utils.h"
#include "framework/physics/physics_namespace.h"
#include "framework/physics/physics_material/multi_group_xs/single_state_mgxs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include <iostream>

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(PhysicsTransportXSCreate, xs, Create);
RegisterLuaFunctionNamespace(PhysicsTransportXSSet, xs, Set);
RegisterLuaFunctionNamespace(PhysicsTransportXSMakeCombined, xs, MakeCombined);
RegisterLuaFunctionNamespace(PhysicsTransportXSSetCombined, xs, SetCombined);
RegisterLuaFunctionNamespace(PhysicsTransportXSGet, xs, Get);
RegisterLuaFunctionNamespace(PhysicsTransportXSExportToOpenSnFormat, xs, ExportToOpenSnFormat);

RegisterLuaConstantAsIs(SINGLE_VALUE, Varying(0));
RegisterLuaConstantAsIs(FROM_ARRAY, Varying(1));
RegisterLuaConstantAsIs(SIMPLEXS0, Varying(20));
RegisterLuaConstantAsIs(SIMPLEXS1, Varying(21));
RegisterLuaConstantAsIs(EXISTING, Varying(22));
RegisterLuaConstantAsIs(OPENSN_XSFILE, Varying(23));

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
  LuaPushTableKey(L, "sigma_t", xs->SigmaTotal());
  LuaPushTableKey(L, "sigma_a", xs->SigmaAbsorption());
  LuaPushTableKey(L, "sigma_f", xs->SigmaFission());
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
PhysicsTransportXSCreate(lua_State* L)
{
  auto xs = std::make_shared<SingleStateMGXS>();
  opensn::multigroup_xs_stack.push_back(xs);

  const size_t index = opensn::multigroup_xs_stack.size() - 1;
  LuaPush(L, index);
  return 1;
}

int
PhysicsTransportXSSet(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  const int num_args = lua_gettop(L);

  if (num_args < 3)
    LuaPostArgAmountError(fname, 3, num_args);

  // Process xs handle
  const auto handle = LuaArg<int>(L, 1);

  // Process operation id
  const auto operation_index = LuaArg<int>(L, 2);

  std::shared_ptr<SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<SingleStateMGXS>(
      opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle));
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
    if (num_args != 4)
      LuaPostArgAmountError(fname, 4, num_args);

    const auto n_grps = LuaArg<int>(L, 3);
    const auto sigma_t = LuaArg<double>(L, 4);

    xs->MakeSimple0(n_grps, sigma_t);
  }
  else if (operation_index == static_cast<int>(OpType::SIMPLEXS1))
  {
    if (num_args != 5)
      LuaPostArgAmountError(fname, 5, num_args);

    const auto n_grps = LuaArg<int>(L, 3);
    const auto sigma_t = LuaArg<double>(L, 4);
    const auto c = LuaArg<double>(L, 5);

    xs->MakeSimple1(n_grps, sigma_t, c);
  }
  else if (operation_index == static_cast<int>(OpType::OPENSN_XSFILE))
  {
    if (num_args != 3)
      LuaPostArgAmountError(fname, 3, num_args);

    auto file_name = LuaArg<std::string>(L, 3);

    xs->MakeFromOpenSnXSFile(file_name);
  }
  else
  {
    opensn::log.LogAllError() << "Unsupported operation in " << fname << ". " << operation_index;
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
}

int
PhysicsTransportXSGet(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  const int num_args = lua_gettop(L);

  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  // Process xs handle
  const auto handle = LuaArg<int>(L, 1);

  std::shared_ptr<SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<SingleStateMGXS>(
      opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle));
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
PhysicsTransportXSMakeCombined(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

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
  auto new_xs = std::make_shared<SingleStateMGXS>();

  new_xs->MakeCombined(combinations);

  opensn::multigroup_xs_stack.push_back(new_xs);
  const lua_Integer num_xs = opensn::multigroup_xs_stack.size();
  LuaPush(L, num_xs - 1);

  return 1;
}

int
PhysicsTransportXSSetCombined(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  int num_args = lua_gettop(L);

  if (num_args < 2)
    LuaPostArgAmountError(fname, 2, num_args);

  // Process xs handle
  const auto xs_handle = LuaArg<int>(L, 1);

  std::shared_ptr<SingleStateMGXS> xs;
  try
  {
    xs = std::dynamic_pointer_cast<SingleStateMGXS>(
      opensn::GetStackItemPtr(opensn::multigroup_xs_stack, xs_handle));
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

  xs->MakeCombined(combinations);

  return 0;
}

int
PhysicsTransportXSExportToOpenSnFormat(lua_State* L)
{
  const auto fname = std::string(__FUNCTION__);
  const int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  // Process xs handle
  const auto handle = LuaArg<int>(L, 1);

  std::shared_ptr<MultiGroupXS> xs;
  try
  {
    xs = opensn::GetStackItemPtr(opensn::multigroup_xs_stack, handle);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "ERROR: Invalid cross section handle in call to " << fname << ".";
    opensn::Exit(EXIT_FAILURE);
  }

  // Process file name
  auto file_name = LuaArg<std::string>(L, 2);

  xs->ExportToOpenSnXSFile(file_name);

  return 0;
}

} // namespace opensnlua
