#pragma once

#include "framework/mesh/chi_mesh.h"
#include "framework/mesh/Cell/cell.h"
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

/**
 * Parent volume mesher class.
 */
class VolumeMesher
{
private:
  MeshContinuumPtr grid_ptr_;
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
  explicit VolumeMesher(VolumeMesherType type);
  virtual ~VolumeMesher() = default;

  /**
   * Sets the grid member of the volume mesher.
   */
  void SetContinuum(MeshContinuumPtr& grid);

  /**
   * Gets the smart-pointer for the grid.
   */
  MeshContinuumPtr& GetContinuum();

  /**
   * Sets grid attributes. This is normally a private member of the grid but this class is a friend.
   */
  void SetGridAttributes(MeshAttributes new_attribs, std::array<size_t, 3> ortho_Nis = {0, 0, 0});

  /**
   * Gets the volume mesher's type.
   */
  VolumeMesherType Type() const;

  /**
   * Obtains the xy partition IDs of a cell.
   * Cell xy_partition ids are obtained from the surface mesher.
   */
  static std::pair<int, int> GetCellXYPartitionID(Cell* cell);

  /**
   * Obtains the xyz partition IDs of a cell.
   * Cell xy_partition ids are obtained from
   * the surface mesher. z id is obtained from the volume mesher.
   */
  static std::tuple<int, int, int> GetCellXYZPartitionID(Cell* cell);

  /**
   * Creates 2D polygon cells for each face of an unpartitioned mesh.
   */
  static void CreatePolygonCells(const UnpartitionedMesh& umesh, MeshContinuumPtr& grid);
  /**
   * Sets material id's using a logical volume.
   */
  static void SetMatIDFromLogical(const LogicalVolume& log_vol, bool sense, int mat_id);

  /**
   * Sets material id's using a logical volume.
   */
  static void
  SetBndryIDFromLogical(const LogicalVolume& log_vol, bool sense, const std::string& bndry_name);

  /**
   * Sets material id's for all cells to the specified material id.
   */
  static void SetMatIDToAll(int mat_id);

#ifdef OPENSN_WITH_LUA
  /**
   * Sets material id's using a lua function. The lua function is called with for each cell with 4
   * arguments, the cell's centroid x,y,z values and the cell's current material id.
   *
   * The lua function's prototype should be:
   * \code
   * function LuaFuncName(x,y,z,id)
   *   --stuff
   * end
   * \endcode
   */
  static void SetMatIDFromLuaFunction(const std::string& lua_fname);

  /**Sets boundary id's using a lua function. The lua function is called for each boundary face
   * with 7 arguments, the face's centroid x,y,z values, the face's normal x,y,z values and the
   * face's current boundary id. The function must return a new_bndry_name (string).
   *
   * The lua function's prototype should be:
   * \code
   * function LuaFuncName(x,y,z,nx,ny,nz,id)
   * --stuff
   * end
   * \endcode
   */
  static void SetBndryIDFromLuaFunction(const std::string& lua_fname);
#endif

  /**
   * Sets boundary numbers on boundaries orthogonal to the cardinal directions as "XMAX", "XMIN",
   * "YMAX", "YMIN", "ZMAX", "ZMIN".
   */
  static void SetupOrthogonalBoundaries();

  virtual void Execute();
};

} // namespace chi_mesh
