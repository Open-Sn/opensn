#include "diffusion_mip.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"

#include "mesh/chi_mesh.h"

#define scdouble static_cast<double>

//###################################################################
/**Still searching for a reference for this.
 *
 * For Polygons:
 * Defined from paper  \n
 * Turcksin B, Ragusa J, "Discontinuous diffusion synthetic acceleration
 * for S_n transport on 2D arbitrary polygonal meshes", Journal of
 * Computational Physics 274, pg 356-369, 2014.\n
 * \n
 * Nv = Number of vertices. If Nv <= 4 then the perimeter parameter
 * should be replaced by edge length.*/
double
lbs::acceleration::DiffusionMIPSolver::HPerpendicular(const chi_mesh::Cell& cell, unsigned int f)
{
  const auto& cell_mapping = sdm_.GetCellMapping(cell);
  double hp;

  const size_t num_faces = cell.faces_.size();
  const size_t num_vertices = cell.vertex_ids_.size();

  const double volume = cell_mapping.CellVolume();
  const double face_area = cell_mapping.FaceArea(f);

  /**Lambda to compute surface area.*/
  auto ComputeSurfaceArea = [&cell_mapping, &num_faces]()
  {
    double surface_area = 0.0;
    for (size_t fr = 0; fr < num_faces; ++fr)
      surface_area += cell_mapping.FaceArea(fr);

    return surface_area;
  };

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell.Type() == chi_mesh::CellType::SLAB) hp = volume / 2.0;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    if (num_faces == 3) hp = 2.0 * volume / face_area;
    else if (num_faces == 4)
      hp = volume / face_area;
    else // Nv > 4
    {
      const double surface_area = ComputeSurfaceArea();

      if (num_faces % 2 == 0) hp = 4.0 * volume / surface_area;
      else
      {
        hp = 2.0 * volume / surface_area;
        hp += sqrt(2.0 * volume / (scdouble(num_faces) * sin(2.0 * M_PI / scdouble(num_faces))));
      }
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    const double surface_area = ComputeSurfaceArea();

    if (num_faces == 4) // Tet
      hp = 3 * volume / surface_area;
    else if (num_faces == 6 && num_vertices == 8) // Hex
      hp = volume / surface_area;
    else // Polyhedron
      hp = 6 * volume / surface_area;
  } // Polyhedron
  else
    throw std::logic_error("lbs::acceleration::DiffusionMIPSolver::HPerpendicular: "
                           "Unsupported cell type in call to HPerpendicular");

  return hp;
}

//###################################################################
/**Maps a face, in a discontinuous sense, using the spatial discretization.*/
int
lbs::acceleration::DiffusionMIPSolver::MapFaceNodeDisc(
  const chi_mesh::Cell& cur_cell,
  const chi_mesh::Cell& adj_cell,
  const std::vector<chi_mesh::Vector3>& cc_node_locs,
  const std::vector<chi_mesh::Vector3>& ac_node_locs,
  size_t ccf,
  size_t acf,
  size_t ccfi,
  double epsilon /*=1.0e-12*/)
{
  const auto& cur_cell_mapping = sdm_.GetCellMapping(cur_cell);
  const auto& adj_cell_mapping = sdm_.GetCellMapping(adj_cell);

  const int i = cur_cell_mapping.MapFaceNode(ccf, ccfi);
  const auto& node_i_loc = cc_node_locs[i];

  const size_t adj_face_num_nodes = adj_cell_mapping.NumFaceNodes(acf);

  for (size_t fj = 0; fj < adj_face_num_nodes; ++fj)
  {
    const int j = adj_cell_mapping.MapFaceNode(acf, fj);
    if ((node_i_loc - ac_node_locs[j]).NormSquare() < epsilon) return j;
  }

  throw std::logic_error(
    "lbs::acceleration::DiffusionMIPSolver::MapFaceNodeDisc: Mapping failure.");
}

//###################################################################
/**Calls a lua function with xyz coordinates.
 * \param L The lua state.
 * \param lua_func_name The name used to define this lua function in the lua
 *                      state.
 * \param xyz The xyz coordinates of the point where the function is called.
 *
 * \return The function evaluation.*/
double
lbs::acceleration::DiffusionMIPSolver::CallLuaXYZFunction(lua_State* L,
                                                          const std::string& lua_func_name,
                                                          const chi_mesh::Vector3& xyz)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::"
                            "CallLuaXYZFunction";
  //============= Load lua function
  lua_getglobal(L, lua_func_name.c_str());

  //============= Error check lua function
  if (not lua_isfunction(L, -1))
    throw std::logic_error(fname + " attempted to access lua-function, " + lua_func_name +
                           ", but it seems the function"
                           " could not be retrieved.");

  //============= Push arguments
  lua_pushnumber(L, xyz.x);
  lua_pushnumber(L, xyz.y);
  lua_pushnumber(L, xyz.z);

  //============= Call lua function
  // 3 arguments, 1 result (double), 0=original error object
  double lua_return;
  if (lua_pcall(L, 3, 1, 0) == 0)
  {
    LuaCheckNumberValue(fname, L, -1);
    lua_return = lua_tonumber(L, -1);
  }
  else
    throw std::logic_error(fname + " attempted to call lua-function, " + lua_func_name +
                           ", but the call failed." + xyz.PrintStr());

  lua_pop(L, 1); // pop the double, or error code

  return lua_return;
}