#pragma once

#include "diffusion.h"
#include "chi_lua.h"

//############################################### Forward declarations
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

//############################################### Namespace lbs::acceleration
namespace lbs::acceleration
{

/**Generalized diffusion solver for both WGDSA and TGDSA based on the MIP-method
 * of Bruno Turcksin and Jean Ragusa.*/
class DiffusionMIPSolver : public lbs::acceleration::DiffusionSolver
{
public:
  // 00
  DiffusionMIPSolver(std::string text_name,
                     const chi_math::SpatialDiscretization& sdm,
                     const chi_math::UnknownManager& uk_man,
                     std::map<uint64_t, BoundaryCondition> bcs,
                     MatID2XSMap map_mat_id_2_xs,
                     const std::vector<UnitCellMatrices>& unit_cell_matrices,
                     bool verbose);

  // 02a
  void AssembleAand_b_wQpoints(const std::vector<double>& q_vector);
  // 02b
  void Assemble_b_wQpoints(const std::vector<double>& q_vector);

  // 02c
  void AssembleAand_b(const std::vector<double>& q_vector) override;
  // 02d
  void Assemble_b(const std::vector<double>& q_vector) override;
  void Assemble_b(Vec petsc_q_vector) override;

  // 05
  double HPerpendicular(const chi_mesh::Cell& cell, unsigned int f);

  int MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                      const chi_mesh::Cell& adj_cell,
                      const std::vector<chi_mesh::Vector3>& cc_node_locs,
                      const std::vector<chi_mesh::Vector3>& ac_node_locs,
                      size_t ccf,
                      size_t acf,
                      size_t ccfi,
                      double epsilon = 1.0e-12);
  static double
  CallLuaXYZFunction(lua_State* L, const std::string& lua_func_name, const chi_mesh::Vector3& xyz);

  virtual ~DiffusionMIPSolver() = default;
};

} // namespace lbs::acceleration
