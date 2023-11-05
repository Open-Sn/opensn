#include "opensn/framework/chi_lua.h"

#include "opensn/framework/mesh/MeshHandler/chi_meshhandler.h"
#include "opensn/framework/mesh/SurfaceMesher/surfacemesher.h"
#include "opensn/framework/mesh/VolumeMesher/Extruder/volmesher_extruder.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/mesh/LogicalVolume/LogicalVolume.h"

#include <iostream>
#include "volumemesher_lua.h"
#include "opensn/framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiVolumeMesherSetProperty);

RegisterLuaConstantAsIs(FORCE_POLYGONS, chi_data_types::Varying(1));
RegisterLuaConstantAsIs(MESH_GLOBAL, chi_data_types::Varying(2));
RegisterLuaConstantAsIs(PARTITION_Z, chi_data_types::Varying(3));
RegisterLuaConstantAsIs(VOLUMEPARTITION_Y, chi_data_types::Varying(4));
RegisterLuaConstantAsIs(VOLUMEPARTITION_X, chi_data_types::Varying(5));
RegisterLuaConstantAsIs(CUTS_Z, chi_data_types::Varying(6));
RegisterLuaConstantAsIs(CUTS_Y, chi_data_types::Varying(7));
RegisterLuaConstantAsIs(CUTS_X, chi_data_types::Varying(8));
RegisterLuaConstantAsIs(PARTITION_TYPE, chi_data_types::Varying(9));
RegisterLuaConstantAsIs(KBA_STYLE_XYZ, chi_data_types::Varying(2));
RegisterLuaConstantAsIs(PARMETIS, chi_data_types::Varying(3));
RegisterLuaConstantAsIs(EXTRUSION_LAYER, chi_data_types::Varying(10));
RegisterLuaConstantAsIs(MATID_FROMLOGICAL, chi_data_types::Varying(11));
RegisterLuaConstantAsIs(BNDRYID_FROMLOGICAL, chi_data_types::Varying(12));
RegisterLuaConstantAsIs(MATID_FROM_LUA_FUNCTION, chi_data_types::Varying(13));
RegisterLuaConstantAsIs(BNDRYID_FROM_LUA_FUNCTION, chi_data_types::Varying(14));

RegisterLuaFunctionAsIs(chiVolumeMesherSetKBAPartitioningPxPyPz);
RegisterLuaFunctionAsIs(chiVolumeMesherSetKBACutsX);
RegisterLuaFunctionAsIs(chiVolumeMesherSetKBACutsY);
RegisterLuaFunctionAsIs(chiVolumeMesherSetKBACutsZ);

// #############################################################################
/** Sets a volume mesher property.

\param PropertyIndex int Index of the property to change. See below
\param PropertyValue varying Value of the property.

##_

###PropertyIndex:
 FORCE_POLYGONS = <B>PropertyValue:[bool]</B> Forces the 2D Meshing to use
                  polygon cells even if the underlying surface mesh is
                  triangles. Expects a boolean value.\n
 MESH_GLOBAL = <B>PropertyValue:[bool]</B> Generate/Read the full mesh at each
               location. Expects a boolean value [Default=true].\n
 PARTITION_TYPE = <B>PartitionType</B>. See below.\n
 EXTRUSION_LAYER = <B>PropertyValue:[double,(int),(char)]</B> Adds a layer to
the extruder volume mesher if it exists. Expects 1 required parameter, the layer
height, followed by 2 optional parameters: number of subdivisions (defaults to
1), and layer id (char)(defaults to nothing). Only supported if partition-type
is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```. \n
 CUTS_X = Adds a cut at the given x-value. Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 CUTS_Y = Adds a cut at the given y-value. Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 CUTS_Z = Adds a cut at the given z-value. Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 PARTITION_X   = <B>PropertyValue:[int]</B> Number of partitions in X.
                    Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 PARTITION_Y   = <B>PropertyValue:[int]</B> Number of partitions in Y.
                    Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 PARTITION_Z   = <B>PropertyValue:[int]</B> Number of partitions in Z.
                    Only supported if partition-type is
                   ```KBA_STYLE_XY``` or ```KBA_STYLE_XYZ```.\n
 MATID_FROMLOGICAL = <B>LogicalVolumeHandle:[int],Mat_id:[int],
                     Sense:[bool](Optional, default:true)</B> Sets the material
                     id of cells that meet the sense requirement for the given
                     logical volume.\n
 BNDRYID_FROMLOGICAL = <B>LogicalVolumeHandle:[int],Bndry_name:[string],
                     Sense:[bool](Optional, default:true)</B> Sets the cell
                     boundary id to the specified value for cells
                     that meet the sense requirement for the given
                     logical volume.\n
 MATID_FROM_LUA_FUNCTION = <B>LuaFunctionName:[string]</B>. For each cell, will
                           call a lua function that can change the material id.
                           The lua function must have 4 parameters,
                           the cell's centroid x,y,z values (doubles) and the
                           current cell-material id (int). The function must
                           return a material id.
 BNDRYID_FROM_LUA_FUNCTION = <B>LuaFunctionName:[string]</B>. For each boundary
                           face, will call a lua function that can change the
                           boundary id.
                           The lua function must have 7 parameters,
                           the face's centroid x,y,z values (doubles), the
                           face's normal x,y,z values (double), and the
                           current face-boundary id (int). The function must
                           return a boundary id.
## _

### PartitionType
Can be any of the following:
 - KBA_STYLE_XYZ
 - PARMETIS

\ingroup LuaVolumeMesher
\author Jan*/
int
chiVolumeMesherSetProperty(lua_State* L)
{
  const std::string fname = "chiVolumeMesherSetProperty";
  //============================================= Get current mesh handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();
  auto& volume_mesher = cur_hndlr.GetVolumeMesher();

  //============================================= Get property index
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  int property_index = lua_tonumber(L, 1);

  typedef chi_mesh::VolumeMesherProperty VMP;

  //============================================= Selects property
  if (property_index == VMP::FORCE_POLYGONS)
  {
    bool force_condition = lua_toboolean(L, 2);
    volume_mesher.options.force_polygons = force_condition;
  }

  else if (property_index == VMP::MESH_GLOBAL)
  {
    bool mesh_global = lua_toboolean(L, 2);
    volume_mesher.options.mesh_global = mesh_global;
  }

  else if (property_index == VMP::PARTITION_Z)
  {
    int pz = lua_tonumber(L, 2);
    volume_mesher.options.partition_z = pz;
    Chi::log.LogAllVerbose1() << "Partition z set to " << pz;
  }
  else if (property_index == VMP::PARTITION_Y)
  {
    int p = lua_tonumber(L, 2);
    volume_mesher.options.partition_y = p;
    Chi::log.LogAllVerbose1() << "Partition y set to " << p;
  }
  else if (property_index == VMP::PARTITION_X)
  {
    int p = lua_tonumber(L, 2);
    volume_mesher.options.partition_x = p;
    Chi::log.LogAllVerbose1() << "Partition x set to " << p;
  }
  else if (property_index == VMP::CUTS_Z)
  {
    double p = lua_tonumber(L, 2);
    volume_mesher.options.zcuts.push_back(p);
  }
  else if (property_index == VMP::CUTS_Y)
  {
    double p = lua_tonumber(L, 2);
    cur_hndlr.GetSurfaceMesher().AddYCut(p);
    volume_mesher.options.ycuts.push_back(p);
  }
  else if (property_index == VMP::CUTS_X)
  {
    double p = lua_tonumber(L, 2);
    cur_hndlr.GetSurfaceMesher().AddXCut(p);
    volume_mesher.options.xcuts.push_back(p);
  }
  else if (property_index == VMP::PARTITION_TYPE)
  {
    int p = lua_tonumber(L, 2);
    if (p >= chi_mesh::VolumeMesher::PartitionType::KBA_STYLE_XYZ and
        p <= chi_mesh::VolumeMesher::PartitionType::PARMETIS)
      volume_mesher.options.partition_type = (chi_mesh::VolumeMesher::PartitionType)p;
    else
    {
      Chi::log.LogAllError() << "Unsupported partition type used in call to " << fname << ".";
      Chi::Exit(EXIT_FAILURE);
    }
  }

  else if (property_index == VMP::EXTRUSION_LAYER)
  {
    if (volume_mesher.Type() == chi_mesh::VolumeMesherType::EXTRUDER)
    {
      auto& mesher = dynamic_cast<chi_mesh::VolumeMesherExtruder&>(volume_mesher);

      double layer_height = lua_tonumber(L, 2);
      int subdivisions = 1;

      if (num_args >= 3) { subdivisions = lua_tonumber(L, 3); }
      chi_mesh::VolumeMesherExtruder::MeshLayer new_layer;
      new_layer.height = layer_height;
      new_layer.sub_divisions = subdivisions;

      if (num_args == 4) { new_layer.name = std::string(lua_tostring(L, 4)); }
      mesher.AddLayer(new_layer);
    }
    else
    {
      fprintf(stdout,
              "VolumeMesherExtruder is not the current volume mesher"
              " therefore the z-layer property is ignored.\n");
    }
  }

  else if (property_index == VMP::MATID_FROMLOGICAL)
  {
    if (!((num_args == 3) || (num_args == 4)))
    {
      Chi::log.LogAllError() << "Invalid amount of arguments used for "
                                "chiVolumeMesherSetProperty("
                                "MATID_FROMLOGICAL...";
      Chi::Exit(EXIT_FAILURE);
    }
    int volume_hndl = lua_tonumber(L, 2);
    int mat_id = lua_tonumber(L, 3);
    int sense = true;
    if (num_args == 4) sense = lua_toboolean(L, 4);

    const auto& log_vol =
      Chi::GetStackItem<chi_mesh::LogicalVolume>(Chi::object_stack, volume_hndl, fname);

    chi_mesh::VolumeMesher::SetMatIDFromLogical(log_vol, sense, mat_id);
  }

  else if (property_index == VMP::BNDRYID_FROMLOGICAL)
  {
    if (!((num_args == 3) || (num_args == 4)))
    {
      Chi::log.LogAllError() << "Invalid amount of arguments used for "
                                "chiVolumeMesherSetProperty("
                                "BNDRYID_FROMLOGICAL...";
      Chi::Exit(EXIT_FAILURE);
    }
    LuaCheckNilValue(fname, L, 2);
    LuaCheckStringValue(fname, L, 3);
    int volume_hndl = lua_tonumber(L, 2);
    std::string bndry_name = lua_tostring(L, 3);
    int sense = true;
    if (num_args == 4) sense = lua_toboolean(L, 4);

    if (bndry_name.empty())
      throw std::invalid_argument(fname + ": argument 3 must not be"
                                          "an empty string.");

    const auto& log_vol =
      Chi::GetStackItem<chi_mesh::LogicalVolume>(Chi::object_stack, volume_hndl, fname);

    chi_mesh::VolumeMesher::SetBndryIDFromLogical(log_vol, sense, bndry_name);
  }
  else if (property_index == VMP::MATID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    chi_mesh::VolumeMesher::SetMatIDFromLuaFunction(lua_fname);
  }
  else if (property_index == VMP::BNDRYID_FROM_LUA_FUNCTION)
  {
    LuaCheckStringValue(fname, L, 2);

    const std::string lua_fname = lua_tostring(L, 2);

    chi_mesh::VolumeMesher::SetBndryIDFromLuaFunction(lua_fname);
  }
  else
  {
    Chi::log.LogAllError() << "Invalid property specified " << property_index
                           << " in call to chiVolumeMesherSetProperty().";
    Chi::Exit(EXIT_FAILURE);
  }

  return 0;
}

// ###################################################################
/**Sets the Px, Py and Pz partititioning parameters for a
 * KBA-type partitioning. This also fixes the process count required to
 * a total of Px*Py*Pz.
 *
\param Px int Number partitions in x.
\param Py int Number partitions in y.
\param Pz int Number partitions in z.
\ingroup LuaVolumeMesher
 */
int
chiVolumeMesherSetKBAPartitioningPxPyPz(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(__FUNCTION__, 3, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);
  LuaCheckNilValue(__FUNCTION__, L, 2);
  LuaCheckNilValue(__FUNCTION__, L, 3);

  //============================================= Get current mesh handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();
  auto& vol_mesher = cur_hndlr.GetVolumeMesher();

  int px = lua_tonumber(L, 1);
  int py = lua_tonumber(L, 2);
  int pz = lua_tonumber(L, 3);

  vol_mesher.options.partition_x = px;
  vol_mesher.options.partition_y = py;
  vol_mesher.options.partition_z = pz;

  return 0;
}

// ###################################################################
/**Sets the x-cuts for KBA type partitioning with a lua array.
\ingroup LuaVolumeMesher
 */
int
chiVolumeMesherSetKBACutsX(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  std::vector<double> cuts;
  LuaPopulateVectorFrom1DArray(__FUNCTION__, L, 1, cuts);

  auto& mesh_handler = chi_mesh::GetCurrentHandler();
  mesh_handler.GetVolumeMesher().options.xcuts = cuts;

  return 0;
}

// ###################################################################
/**Sets the y-cuts for KBA type partitioning with a lua array.
\ingroup LuaVolumeMesher
 */
int
chiVolumeMesherSetKBACutsY(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  std::vector<double> cuts;
  LuaPopulateVectorFrom1DArray(__FUNCTION__, L, 1, cuts);

  auto& mesh_handler = chi_mesh::GetCurrentHandler();
  mesh_handler.GetVolumeMesher().options.ycuts = cuts;

  return 0;
}

// ###################################################################
/**Sets the z-cuts for KBA type partitioning with a lua array.
\ingroup LuaVolumeMesher
 */
int
chiVolumeMesherSetKBACutsZ(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  std::vector<double> cuts;
  LuaPopulateVectorFrom1DArray(__FUNCTION__, L, 1, cuts);

  auto& mesh_handler = chi_mesh::GetCurrentHandler();
  mesh_handler.GetVolumeMesher().options.zcuts = cuts;

  return 0;
}
