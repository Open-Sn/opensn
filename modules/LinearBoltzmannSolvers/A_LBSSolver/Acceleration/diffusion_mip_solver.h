#pragma once

#include "opensn/modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/diffusion.h"
#include "opensn/framework/chi_lua.h"

// ############################################### Forward declarations
namespace chi_mesh
{
class MeshContinuum;
class Cell;
struct Vector3;
} // namespace chi_mesh

namespace chi_math
{
class SpatialDiscretization;
}

namespace lbs
{
struct UnitCellMatrices;
}

// ############################################### Namespace lbs::acceleration
namespace lbs::acceleration
{

/**
 * Generalized diffusion solver for both WGDSA and TGDSA based on the MIP-method
 * of Bruno Turcksin and Jean Ragusa.
 */
class DiffusionMIPSolver : public lbs::acceleration::DiffusionSolver
{
public:
  DiffusionMIPSolver(std::string text_name,
                     const chi_math::SpatialDiscretization& sdm,
                     const chi_math::UnknownManager& uk_man,
                     std::map<uint64_t, BoundaryCondition> bcs,
                     MatID2XSMap map_mat_id_2_xs,
                     const std::vector<UnitCellMatrices>& unit_cell_matrices,
                     bool verbose);
  virtual ~DiffusionMIPSolver() = default;

  /**
   * Assembles both the matrix and the RHS using quadrature points. These routines exist for
   * implementing MMS.
   */
  void AssembleAand_b_wQpoints(const std::vector<double>& q_vector);

  /**
   * Assembles just the RHS using quadrature points. These routines exist for implementing MMS.
   */
  void Assemble_b_wQpoints(const std::vector<double>& q_vector);

  /**
   * Assembles both the matrix and the RHS using unit cell-matrices. These are the routines used in
   * the production versions.
   */
  void AssembleAand_b(const std::vector<double>& q_vector) override;

  /**
   * Assembles the RHS using unit cell-matrices. These are the routines used in the production
   * versions.
   */
  void Assemble_b(const std::vector<double>& q_vector) override;
  void Assemble_b(Vec petsc_q_vector) override;

  /**
   * Still searching for a reference for this.
   *
   * For Polygons:
   * Defined from paper  \n
   * Turcksin B, Ragusa J, "Discontinuous diffusion synthetic acceleration
   * for S_n transport on 2D arbitrary polygonal meshes", Journal of
   * Computational Physics 274, pg 356-369, 2014.\n
   * \n
   * Nv = Number of vertices. If Nv <= 4 then the perimeter parameter
   * should be replaced by edge length.
   */
  double HPerpendicular(const chi_mesh::Cell& cell, unsigned int f);

  /**
   * Maps a face, in a discontinuous sense, using the spatial discretization.
   */
  int MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                      const chi_mesh::Cell& adj_cell,
                      const std::vector<chi_mesh::Vector3>& cc_node_locs,
                      const std::vector<chi_mesh::Vector3>& ac_node_locs,
                      size_t ccf,
                      size_t acf,
                      size_t ccfi,
                      double epsilon = 1.0e-12);

  /**
   * Calls a lua function with xyz coordinates.
   * \param L The lua state.
   * \param lua_func_name The name used to define this lua function in the lua
   *                      state.
   * \param xyz The xyz coordinates of the point where the function is called.
   *
   * \return The function evaluation.
   */
  static double
  CallLuaXYZFunction(lua_State* L, const std::string& lua_func_name, const chi_mesh::Vector3& xyz);
};

} // namespace lbs::acceleration
