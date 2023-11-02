#include "chi_modules_lua.h"

#include "DiffusionSolver/lua/diffusion_lua.h"
#include "CFEMDiffusion/lua/ds_lua_utils.h"
#include "DFEMDiffusion/lua/ip_lua_utils.h"
#include "FVDiffusion/lua/ds_lua_utils.h"
#include "MGDiffusion/lua/mgds_lua_utils.h"
#include "A_LBSSolver/lua/lbs_lua_utils.h"

void chi_modules::lua_utils::RegisterLuaEntities(lua_State *L)
{
  lbs::common_lua_utils::RegisterLuaEntities(L);
  diffusion_solver::lua_utils::RegisterLuaEntities(L);

  cfem_diffusion::cfem_diffusion_lua_utils::RegisterLuaEntities(L);
  dfem_diffusion::dfem_diffusion_lua_utils::RegisterLuaEntities(L);

  mg_diffusion::mgd_lua_utils::RegisterLuaEntities(L);
  fv_diffusion::fv_diffusion_lua_utils::RegisterLuaEntities(L);
}