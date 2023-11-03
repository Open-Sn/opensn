#pragma once

#include "mesh/chi_mesh.h"
#include "mesh/Cell/cell.h"
#include <array>

#include <array>

#include <array>

namespace chi_mesh
{
enum class VolumeMesherType
{
  EXTRUDER = 4,
  UNPARTITIONED = 6
};
enum VolumeMesherProperty
{
  FORCE_POLYGONS = 1,
  MESH_GLOBAL = 2,
  PARTITION_Z = 3,
  PARTITION_Y = 4,
  PARTITION_X = 5,
  CUTS_Z = 6,
  CUTS_Y = 7,
  CUTS_X = 8,
  PARTITION_TYPE = 9,
  EXTRUSION_LAYER = 10,
  MATID_FROMLOGICAL = 11,
  BNDRYID_FROMLOGICAL = 12,
  MATID_FROM_LUA_FUNCTION = 13,
  BNDRYID_FROM_LUA_FUNCTION = 14
};
} // namespace chi_mesh

//######################################################### Class def
/**Parent volume mesher class.*/
class chi_mesh::VolumeMesher
{
private:
  chi_mesh::MeshContinuumPtr grid_ptr_;
  const VolumeMesherType type_;

public:
  enum PartitionType
  {
    KBA_STYLE_XYZ = 2,
    PARMETIS = 3
  };
  struct VOLUME_MESHER_OPTIONS
  {
    bool force_polygons = true; // TODO: Remove this option
    bool mesh_global = false;   // TODO: Remove this option
    int partition_x = 1;
    int partition_y = 1;
    int partition_z = 1;

    std::vector<double> xcuts;
    std::vector<double> ycuts;
    std::vector<double> zcuts;
    PartitionType partition_type = PARMETIS;
  };
  VOLUME_MESHER_OPTIONS options;

public:
  explicit
    // 01 Utils
    VolumeMesher(VolumeMesherType type);
  virtual ~VolumeMesher() = default;
  void SetContinuum(MeshContinuumPtr& grid);
  MeshContinuumPtr& GetContinuum();
  void SetGridAttributes(MeshAttributes new_attribs, std::array<size_t, 3> ortho_Nis = {0, 0, 0});
  VolumeMesherType Type() const;

  // 01a
  static std::pair<int, int> GetCellXYPartitionID(chi_mesh::Cell* cell);
  static std::tuple<int, int, int> GetCellXYZPartitionID(chi_mesh::Cell* cell);
  // 01b
  static void CreatePolygonCells(const chi_mesh::UnpartitionedMesh& umesh,
                                 chi_mesh::MeshContinuumPtr& grid);
  // 01c
  static void SetMatIDFromLogical(const chi_mesh::LogicalVolume& log_vol, bool sense, int mat_id);
  static void SetBndryIDFromLogical(const chi_mesh::LogicalVolume& log_vol,
                                    bool sense,
                                    const std::string& bndry_name);
  static void SetMatIDToAll(int mat_id);

  static void SetMatIDFromLuaFunction(const std::string& lua_fname);
  static void SetBndryIDFromLuaFunction(const std::string& lua_fname);
  // 01d
  static void SetupOrthogonalBoundaries();
  // 02
  virtual void Execute();
};

