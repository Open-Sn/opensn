// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "framework/data_types/varying.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "mesh_ortho_macros.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaConstantInNamespace(OrthoBoundaryID, XMIN, Varying(static_cast<int>(XMIN)));
RegisterLuaConstantInNamespace(OrthoBoundaryID, XMAX, Varying(static_cast<int>(XMAX)));
RegisterLuaConstantInNamespace(OrthoBoundaryID, YMIN, Varying(static_cast<int>(YMIN)));
RegisterLuaConstantInNamespace(OrthoBoundaryID, YMAX, Varying(static_cast<int>(YMAX)));
RegisterLuaConstantInNamespace(OrthoBoundaryID, ZMIN, Varying(static_cast<int>(ZMIN)));
RegisterLuaConstantInNamespace(OrthoBoundaryID, ZMAX, Varying(static_cast<int>(ZMAX)));

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
  {
    for (auto& face : cell.faces_)
    {
      if (not face.has_neighbor_)
      {
        Vector3& n = face.normal_;

        std::string boundary_name;
        if (n.Dot(ihat) < -0.999)
          boundary_name = "XMIN";
        else if (n.Dot(ihat) > 0.999)
          boundary_name = "XMAX";
        else if (n.Dot(jhat) < -0.999)
          boundary_name = "YMIN";
        else if (n.Dot(jhat) > 0.999)
          boundary_name = "YMAX";
        else if (n.Dot(khat) < -0.999)
          boundary_name = "ZMIN";
        else if (n.Dot(khat) > 0.999)
          boundary_name = "ZMAX";

        uint64_t bndry_id = vol_cont->MakeBoundaryID(boundary_name);

        face.neighbor_id_ = bndry_id;

        vol_cont->GetBoundaryIDMap()[bndry_id] = boundary_name;
      }
    }
  }

  opensn::mpi_comm.barrier();
  opensn::log.Log() << program_timer.GetTimeString() << " Done setting orthogonal boundaries.";

  return LuaReturn(L);
}

} // namespace opensnlua
