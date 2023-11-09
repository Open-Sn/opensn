#include "framework/lua.h"

#include "framework/mesh/logical_volume/logical_volume.h"

#include "framework/object_factory.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/runtime.h"

#include "framework/logging/log.h"

#define scint static_cast<int>

#include "logic_volume_lua.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiLogicalVolumeCreate);
RegisterLuaConstantAsIs(SPHERE, chi_data_types::Varying(1));
RegisterLuaConstantAsIs(SPHERE_ORIGIN, chi_data_types::Varying(2));
RegisterLuaConstantAsIs(RPP, chi_data_types::Varying(3));
RegisterLuaConstantAsIs(RCC, chi_data_types::Varying(4));
RegisterLuaConstantAsIs(SURFACE, chi_data_types::Varying(9));
RegisterLuaConstantAsIs(BOOLEAN, chi_data_types::Varying(10));
RegisterLuaFunctionAsIs(chiLogicalVolumePointSense);

int
chiLogicalVolumeCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  const int num_args = lua_gettop(L);
  const int type_index = lua_tonumber(L, 1);

  enum class LogicalVolumeType
  {
    LVSPHERE = 1,
    LVSPHERE_ORIGIN = 2,
    LVRPP = 3,
    LVRCC = 4,
    LVSURFACE = 9,
    LVBOOLEAN = 10
  };

  const int LVSPHERE = scint(LogicalVolumeType::LVSPHERE);
  const int LVSPHERE_ORIGIN = scint(LogicalVolumeType::LVSPHERE_ORIGIN);
  const int LVRPP = scint(LogicalVolumeType::LVRPP);
  const int LVRCC = scint(LogicalVolumeType::LVRCC);
  const int LVSURFACE = scint(LogicalVolumeType::LVSURFACE);
  const int LVBOOLEAN = scint(LogicalVolumeType::LVBOOLEAN);

  auto& object_maker = ChiObjectFactory::GetInstance();

  // Sphere at origin
  if (type_index == LVSPHERE_ORIGIN)
  {
    if (num_args != 2)
    {
      Chi::log.Log0Error() << "Incorrect amount of arguments provided "
                              "for chiMeshCreateLogicalVolume(SO...";
      Chi::Exit(EXIT_FAILURE);
    }
    double r = lua_tonumber(L, 2);
    // auto log_vol = std::make_shared<chi_mesh::SphereLogicalVolume>(r);
    //
    // chi::object_stack.push_back(log_vol);
    // const size_t index = chi::object_stack.size() - 1;
    // lua_pushinteger(L, static_cast<lua_Integer>(index));

    chi::ParameterBlock params;
    params.AddParameter("r", r);

    const size_t handle =
      object_maker.MakeRegisteredObjectOfType("chi_mesh::SphereLogicalVolume", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }

  // Sphere at arb loc
  else if (type_index == LVSPHERE)
  {
    if (num_args != 5)
    {
      Chi::log.Log0Error() << "Incorrect amount of arguments provided "
                              "for chiMeshCreateLogicalVolume(S...";
      Chi::Exit(EXIT_FAILURE);
    }
    double x = lua_tonumber(L, 2);
    double y = lua_tonumber(L, 3);
    double z = lua_tonumber(L, 4);
    double r = lua_tonumber(L, 5);
    // auto log_vol = std::make_shared<chi_mesh::SphereLogicalVolume>(x, y, z,
    // r);
    //
    // chi::object_stack.push_back(log_vol);
    // const size_t index = chi::object_stack.size() - 1;
    // lua_pushinteger(L, static_cast<lua_Integer>(index));

    chi::ParameterBlock params;
    params.AddParameter("r", r);
    params.AddParameter("x", x);
    params.AddParameter("y", y);
    params.AddParameter("z", z);

    const size_t handle =
      object_maker.MakeRegisteredObjectOfType("chi_mesh::SphereLogicalVolume", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }

  // RPP
  else if (type_index == LVRPP)
  {
    if (num_args != 7)
    {
      Chi::log.Log0Error() << "Incorrect amount of arguments provided "
                              "for chiMeshCreateLogicalVolume(RPP...";
      Chi::Exit(EXIT_FAILURE);
    }
    double xmin = lua_tonumber(L, 2);
    double xmax = lua_tonumber(L, 3);
    double ymin = lua_tonumber(L, 4);
    double ymax = lua_tonumber(L, 5);
    double zmin = lua_tonumber(L, 6);
    double zmax = lua_tonumber(L, 7);
    // auto log_vol = std::make_shared<chi_mesh::RPPLogicalVolume>(
    //   xmin, xmax, ymin, ymax, zmin, zmax);
    //
    // chi::object_stack.push_back(log_vol);
    // const size_t index = chi::object_stack.size() - 1;
    // lua_pushinteger(L, static_cast<lua_Integer>(index));

    chi::ParameterBlock params;
    params.AddParameter("xmin", xmin);
    params.AddParameter("xmax", xmax);
    params.AddParameter("ymin", ymin);
    params.AddParameter("ymax", ymax);
    params.AddParameter("zmin", zmin);
    params.AddParameter("zmax", zmax);

    const size_t handle =
      object_maker.MakeRegisteredObjectOfType("chi_mesh::RPPLogicalVolume", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }

  // RCC
  else if (type_index == LVRCC)
  {
    if (num_args != 8)
    {
      Chi::log.Log0Error() << "Incorrect amount of arguments provided "
                              "for chiMeshCreateLogicalVolume(RCC...";
      Chi::Exit(EXIT_FAILURE);
    }
    double x0 = lua_tonumber(L, 2);
    double y0 = lua_tonumber(L, 3);
    double z0 = lua_tonumber(L, 4);
    double vx = lua_tonumber(L, 5);
    double vy = lua_tonumber(L, 6);
    double vz = lua_tonumber(L, 7);
    double r = lua_tonumber(L, 8);
    // auto log_vol =
    //   std::make_shared<chi_mesh::RCCLogicalVolume>(x0, y0, z0, vx, vy, vz,
    //   r);
    //
    // chi::object_stack.push_back(log_vol);
    // const size_t index = chi::object_stack.size() - 1;
    // lua_pushinteger(L, static_cast<lua_Integer>(index));

    chi::ParameterBlock params;
    params.AddParameter("x0", x0);
    params.AddParameter("y0", y0);
    params.AddParameter("z0", z0);
    params.AddParameter("vx", vx);
    params.AddParameter("vy", vy);
    params.AddParameter("vz", vz);
    params.AddParameter("r", r);

    const size_t handle =
      object_maker.MakeRegisteredObjectOfType("chi_mesh::RCCLogicalVolume", params);

    Chi::log.Log0Verbose1() << "Created RCC Logical volume with x0,y0,z0,vx,vy,vz,r = " << x0 << " "
                            << y0 << " " << z0 << " " << vx << " " << vy << " " << vz << " " << r;

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }
  // SURFACE
  else if (type_index == LVSURFACE)
  {
    if (num_args != 2) LuaPostArgAmountError("chiMeshCreateLogicalVolume:SURFACE", 2, num_args);

    int surf_mesh_hndle = lua_tonumber(L, 2);

    // auto surf_mesh_ptr = chi::GetStackItemPtr<chi_mesh::SurfaceMesh>(
    //   chi::surface_mesh_stack, surf_mesh_hndle, fname);
    //
    // auto log_vol =
    //   std::make_shared<chi_mesh::SurfaceMeshLogicalVolume>(surf_mesh_ptr);
    //
    // chi::object_stack.push_back(log_vol);
    // const size_t index = chi::object_stack.size() - 1;
    // lua_pushinteger(L, static_cast<lua_Integer>(index));
    chi::ParameterBlock params;
    params.AddParameter("surface_mesh_handle", surf_mesh_hndle);

    const size_t handle =
      object_maker.MakeRegisteredObjectOfType("chi_mesh::SurfaceMeshLogicalVolume", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }
  // BOOLEAN
  else if (type_index == LVBOOLEAN)
  {
    // if (num_args % 2 != 0)
    //{
    //   chi::log.Log0Error() << "Incorrect amount of arguments provided "
    //                           "for chiMeshCreateLogicalVolume(BOOLEAN..."
    //                           " Expected pairs of (bool,volumeHandle)";
    //   chi::Exit(EXIT_FAILURE);
    // }

    // auto bool_vol = std::make_shared<chi_mesh::BooleanLogicalVolume>();
    //
    // int num_pairs = num_args / 2;
    // for (int p = 0; p < num_pairs; p++)
    //{
    //   // Checking first part of pair
    //   if (not lua_isboolean(L, 2 * p))
    //   {
    //     chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
    //                             " argument "
    //                          << 2 * p
    //                          << " expected to be "
    //                             "Boolean. Found not to be";
    //     chi::Exit(EXIT_FAILURE);
    //   }
    //   // Checking second part of pair
    //   if (not lua_isnumber(L, 2 * p + 1))
    //   {
    //     chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
    //                             " argument "
    //                          << 2 * p + 1
    //                          << " expected to be "
    //                             "number. Found not to be";
    //     chi::Exit(EXIT_FAILURE);
    //   }
    //   if (lua_tonumber(L, 2 * p + 1) >=
    //       static_cast<lua_Number>(chi::object_stack.size()))
    //   {
    //     chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
    //                             " argument "
    //                          << 2 * p + 1 << " points to non-existent
    //                          volume.";
    //     chi::Exit(EXIT_FAILURE);
    //   }
    //
    //   bool logic = lua_toboolean(L, 2 * p);
    //   bool handle = lua_tonumber(L, 2 * p + 1);
    //
    //   auto obj_ptr = chi::GetStackItemPtr(chi::object_stack, handle, fname);
    //
    //   typedef std::shared_ptr<chi_mesh::LogicalVolume> LogVolPtr;
    //   auto p_ref_vol =
    //     std::dynamic_pointer_cast<chi_mesh::LogicalVolume>(obj_ptr);
    //
    //   std::pair<bool, LogVolPtr> combo(logic, p_ref_vol);
    //
    //   bool_vol->parts.push_back(combo);
    // }
    //
    // chi::object_stack.push_back(bool_vol);
    // const size_t index = chi::object_stack.size() - 1;
    // lua_pushinteger(L, static_cast<lua_Integer>(index));

    ChiInvalidArgumentIf(num_args % 2 != 0,
                         "Incorrect amount of arguments provided for "
                         "chiMeshCreateLogicalVolume(BOOLEAN..."
                         " Expected pairs of (bool,volumeHandle)");

    chi::ParameterBlock params;

    int num_pairs = num_args / 2;
    for (int p = 0; p < num_pairs; p++)
    {
      // Checking first part of pair
      if (not lua_isboolean(L, 2 * p))
      {
        Chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                " argument "
                             << 2 * p
                             << " expected to be "
                                "Boolean. Found not to be";
        Chi::Exit(EXIT_FAILURE);
      }
      // Checking second part of pair
      if (not lua_isnumber(L, 2 * p + 1))
      {
        Chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                " argument "
                             << 2 * p + 1
                             << " expected to be "
                                "number. Found not to be";
        Chi::Exit(EXIT_FAILURE);
      }
      if (lua_tointeger(L, 2 * p + 1) >= static_cast<lua_Number>(Chi::object_stack.size()))
      {
        Chi::log.Log0Error() << "chiMeshCreateLogicalVolume(BOOLEAN..."
                                " argument "
                             << 2 * p + 1 << " points to non-existent volume.";
        Chi::Exit(EXIT_FAILURE);
      }

      const bool logic = lua_toboolean(L, 2 * p);
      const bool handle = lua_tonumber(L, 2 * p + 1);

      chi::ParameterBlock block;
      block.AddParameter("op", logic);
      block.AddParameter("lv", handle);

      params.AddParameter(block);
    }
    params.ChangeToArray();

    const size_t handle =
      object_maker.MakeRegisteredObjectOfType("chi_mesh::BooleanLogicalVolume", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }

  // Unrecognized option
  else
  {
    Chi::log.Log0Error() << "Unrecognized volume type used in "
                            "chiLogicalVolumeCreate.";
    Chi::Exit(EXIT_FAILURE);
  }

  return 1;
}

int
chiLogicalVolumePointSense(lua_State* L)
{
  const std::string fname = "chiLogicalVolumePointSense";
  const int num_args = lua_gettop(L);
  if (num_args != 4) LuaPostArgAmountError(fname, 4, num_args);

  LuaCheckNilValue(fname, L, 1);

  const int lv_handle = lua_tointeger(L, 1);

  const auto& lv = Chi::GetStackItem<chi_mesh::LogicalVolume>(Chi::object_stack, lv_handle, fname);

  const chi_mesh::Vector3 point(lua_tonumber(L, 2), lua_tonumber(L, 3), lua_tonumber(L, 4));

  if (lv.Inside(point)) lua_pushboolean(L, true);
  else
    lua_pushboolean(L, false);

  return 1;
}
