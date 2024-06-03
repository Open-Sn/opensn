// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "framework/data_types/varying.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "lua/framework/mesh/ortho_grids/mesh_ortho_macros.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaConstantInNamespace(OrthoBoundaryID, XMAX, Varying(0));
RegisterLuaConstantInNamespace(OrthoBoundaryID, XMIN, Varying(1));
RegisterLuaConstantInNamespace(OrthoBoundaryID, YMAX, Varying(2));
RegisterLuaConstantInNamespace(OrthoBoundaryID, YMIN, Varying(3));
RegisterLuaConstantInNamespace(OrthoBoundaryID, ZMAX, Varying(4));
RegisterLuaConstantInNamespace(OrthoBoundaryID, ZMIN, Varying(5));

RegisterLuaFunctionInNamespace(MeshSetupOrthogonalBoundaries, mesh, SetupOrthogonalBoundaries);

int
MeshSetupOrthogonalBoundaries(lua_State* L)
{
  opensn::log.Log() << program_timer.GetTimeString() << " Setting orthogonal boundaries.";

  auto vol_cont = GetCurrentMesh();

  const Vector3 ihat(1.0, 0.0, 0.0);
  const Vector3 jhat(0.0, 1.0, 0.0);
  const Vector3 khat(0.0, 0.0, 1.0);

  for (auto& cell : vol_cont->local_cells)
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        Vector3& n = face.normal_;

        std::string boundary_name;
        if (n.Dot(ihat) > 0.999)
          boundary_name = "XMAX";
        else if (n.Dot(ihat) < -0.999)
          boundary_name = "XMIN";
        else if (n.Dot(jhat) > 0.999)
          boundary_name = "YMAX";
        else if (n.Dot(jhat) < -0.999)
          boundary_name = "YMIN";
        else if (n.Dot(khat) > 0.999)
          boundary_name = "ZMAX";
        else if (n.Dot(khat) < -0.999)
          boundary_name = "ZMIN";

        uint64_t bndry_id = vol_cont->MakeBoundaryID(boundary_name);

        face.neighbor_id_ = bndry_id;

        vol_cont->GetBoundaryIDMap()[bndry_id] = boundary_name;
      } // if bndry

  opensn::mpi_comm.barrier();
  opensn::log.Log() << program_timer.GetTimeString() << " Done setting orthogonal boundaries.";

  return LuaReturn(L);
}

} // namespace opensnlua
