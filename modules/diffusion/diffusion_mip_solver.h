// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/diffusion/diffusion.h"
#include "framework/math/functions/function.h"
#include "framework/data_types/vector3.h"
#include <utility>
#include <vector>

namespace opensn
{
class Mesh;
class Cell;
struct Vector3;
class SpatialDiscretization;
struct UnitCellMatrices;

/**
 * Generalized diffusion solver for both WGDSA and TGDSA based on the MIP-method
 * of Bruno Turcksin and Jean Ragusa.
 */
class DiffusionMIPSolver : public DiffusionSolver
{
public:
  DiffusionMIPSolver(std::string name,
                     const opensn::SpatialDiscretization& sdm,
                     const UnknownManager& uk_man,
                     std::map<uint64_t, BoundaryCondition> bcs,
                     MatID2XSMap map_mat_id_2_xs,
                     const std::vector<UnitCellMatrices>& unit_cell_matrices,
                     bool suppress_bcs,
                     bool verbose);
  ~DiffusionMIPSolver() override = default;

  void SetSourceFunction(const ScalarSpatialFunction& function) { source_function_ = function; }

  void SetReferenceSolutionFunction(const ScalarSpatialFunction& function)
  {
    ref_solution_function_ = function;
  }

  /**
   * Sets the angular quadrature (unit directions and weights) of the transport sweep. It is
   * required by the CLOSURE_R boundary type, which uses the quadrature's half-range moment beta_f
   * of each boundary face, and it replaces the constant 1/4 floor of the interior penalty by beta_f
   * of each interior face (the transport operator's own jump dissipation). Call before assembly.
   */
  void SetBoundaryClosureQuadrature(const std::vector<Vector3>& omegas,
                                    const std::vector<double>& weights)
  {
    closure_omegas_ = omegas;
    closure_weights_ = weights;
  }

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
  double HPerpendicular(const Cell& cell, unsigned int f);

  /**
   * Extra factor on the interior-face penalty of polyhedral cells: 2 for the "cell" h_perp
   * convention (OpenSn's original) and for "face_x2", 1 for "face".
   */
  double PolyhedronInteriorPenaltyScale() const;

  /// Maps a face, in a discontinuous sense, using the spatial discretization.
  int MapFaceNodeDisc(const Cell& cur_cell,
                      const Cell& adj_cell,
                      const std::vector<Vector3>& cc_node_locs,
                      const std::vector<Vector3>& ac_node_locs,
                      size_t ccf,
                      size_t acf,
                      size_t ccfi,
                      double epsilon = 1.0e-12);

private:
  /// CLOSURE_R data of a boundary face: half-range moment, orthogonal length, trace constant.
  struct ClosureFaceData
  {
    double beta = 0.0;
    double h = 0.0;
    double c_grad = 0.0;
  };

  /**
   * CLOSURE_R data of every face of `cell` (filled for CLOSURE_R boundary faces only) and the
   * number of such faces.
   */
  std::pair<std::vector<ClosureFaceData>, unsigned int>
  ComputeClosureFaceData(const Cell& cell) const;

  /// Floor of the interior penalty on a face with unit normal `normal`: beta_f or 1/4.
  double PenaltyFloor(const Vector3& normal) const;

  ScalarSpatialFunction source_function_;
  ScalarSpatialFunction ref_solution_function_;
  std::vector<Vector3> closure_omegas_;
  std::vector<double> closure_weights_;
};

} // namespace opensn
