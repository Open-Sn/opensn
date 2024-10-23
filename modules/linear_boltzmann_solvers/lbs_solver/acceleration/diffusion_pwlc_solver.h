// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion.h"

namespace opensn
{

class DiffusionPWLCSolver : public DiffusionSolver
{
public:
  DiffusionPWLCSolver(std::string name,
                      const opensn::SpatialDiscretization& sdm,
                      const UnknownManager& uk_man,
                      std::map<uint64_t, BoundaryCondition> bcs,
                      MatID2XSMap map_mat_id_2_xs,
                      const std::vector<UnitCellMatrices>& unit_cell_matrices,
                      bool suppress_bcs,
                      bool verbose);

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

  /**
   * Assembles the RHS using unit cell-matrices. These are the routines used in the production
   * versions.
   */
  void Assemble_b(Vec petsc_q_vector) override;
};

} // namespace opensn
