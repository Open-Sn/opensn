#include "chi_modules_lua.h"

#include "diffusion_solver/diffusion_lua.h"
#include "cfem_diffusion/ds_lua_utils.h"
#include "dfem_diffusion/ip_lua_utils.h"
#include "fv_diffusion/ds_lua_utils.h"
#include "mg_diffusion/mgds_lua_utils.h"
#include "linear_bolzmann_solvers/lbs_solver/lbs_lua_utils.h"

void
chi_modules::lua_utils::RegisterLuaEntities(lua_State* L)
{
  lbs::common_lua_utils::RegisterLuaEntities(L);
  diffusion_solver::lua_utils::RegisterLuaEntities(L);

  cfem_diffusion::cfem_diffusion_lua_utils::RegisterLuaEntities(L);
  dfem_diffusion::dfem_diffusion_lua_utils::RegisterLuaEntities(L);

  mg_diffusion::mgd_lua_utils::RegisterLuaEntities(L);
  fv_diffusion::fv_diffusion_lua_utils::RegisterLuaEntities(L);
}
