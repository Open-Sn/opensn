// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "gmock/gmock.h"
#include "test/unit/common/mesh_builders.h"

#include "modules/diffusion/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/acceleration.h"

#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

// Homogeneous slab with a uniform source q. Without boundary terms, the MIP form of
// -D phi'' + sigma_r phi = q is solved exactly by the constant phi = q / sigma_r: a constant has
// no gradient and no jump across interior faces. A boundary term moves the solution away from
// that constant near the faces (by 0.87 and 0.54 times it for the two conditions below).
constexpr double slab_length = 4.0;
constexpr unsigned int num_cells = 40;
constexpr double diffusion_coeff = 0.5;
constexpr double sigma_r = 0.25;
constexpr double source = 1.0;
constexpr double phi_no_bc = source / sigma_r;

enum class Assembly
{
  UNIT_CELL_MATRICES, // AssembleAand_b, then Assemble_b
  QUADRATURE_POINTS   // AssembleAand_b_wQpoints, then Assemble_b_wQpoints
};

/// Largest |phi - q / sigma_r| over all DOFs, after the full assembly and after the RHS-only one.
std::pair<double, double>
MaxDeviationFromNoBCSolution(const BoundaryCondition& bc, bool suppress_bcs, Assembly assembly)
{
  auto grid = BuildLineMesh(slab_length, num_cells, 0.0);
  grid->SetUniformBlockID(0);
  grid->SetOrthogonalBoundaries();

  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
  const auto& sdm = *sdm_ptr;
  const auto& one_dof_per_node = sdm.UNITARY_UNKNOWN_MANAGER;
  const size_t num_local_dofs = sdm.GetNumLocalDOFs(one_dof_per_node);

  std::map<uint64_t, BoundaryCondition> bcs;
  for (const auto& [bid, name] : grid->GetBoundaryIDMap())
    bcs[bid] = bc;

  std::map<unsigned int, Multigroup_D_and_sigR> matid_2_xs_map;
  matid_2_xs_map.insert(std::make_pair(0, Multigroup_D_and_sigR{{diffusion_coeff}, {sigma_r}}));

  std::vector<UnitCellMatrices> unit_cell_matrices(grid->GetLocalCellCount());
  for (const auto& cell : grid->GetLocalCells())
    unit_cell_matrices[cell->local_id] = ComputeUnitCellIntegrals(sdm, *cell);

  DiffusionMIPSolver solver("MIPSuppressBCs",
                            sdm,
                            one_dof_per_node,
                            bcs,
                            matid_2_xs_map,
                            unit_cell_matrices,
                            suppress_bcs,
                            false);
  solver.options.residual_tolerance = 1.0e-12;
  solver.options.max_iters = 1000;
  solver.options.perform_symmetry_check = true;
  solver.Initialize();

  const std::vector<double> q_vector(num_local_dofs, source);
  std::vector<double> phi(num_local_dofs, 0.0);

  auto max_deviation = [&phi]()
  {
    double local_max = 0.0;
    for (const double value : phi)
      local_max = std::max(local_max, std::fabs(value - phi_no_bc));
    double global_max = 0.0;
    opensn::mpi_comm.all_reduce(local_max, global_max, mpi::op::max<double>());
    return global_max;
  };

  if (assembly == Assembly::UNIT_CELL_MATRICES)
    solver.AssembleAand_b(q_vector);
  else
    solver.AssembleAand_b_wQpoints(q_vector);
  solver.Solve(phi);
  const double full_assembly_deviation = max_deviation();

  if (assembly == Assembly::UNIT_CELL_MATRICES)
    solver.Assemble_b(q_vector);
  else
    solver.Assemble_b_wQpoints(q_vector);
  solver.Solve(phi);
  const double rhs_assembly_deviation = max_deviation();

  return {full_assembly_deviation, rhs_assembly_deviation};
}

void
CheckSuppressedBCs(Assembly assembly)
{
  // Both conditions add matrix and right-hand-side terms.
  const std::vector<std::pair<std::string, BoundaryCondition>> conditions = {
    {"Dirichlet", {BCType::DIRICHLET, {0.5, 0.0, 0.0}}},
    {"Robin", {BCType::ROBIN, {0.25, 0.5, 0.125}}}};

  for (const auto& [label, bc] : conditions)
  {
    SCOPED_TRACE(label);

    const auto [applied_full, applied_rhs] = MaxDeviationFromNoBCSolution(bc, false, assembly);
    EXPECT_GT(applied_full, 0.1 * phi_no_bc);
    EXPECT_GT(applied_rhs, 0.1 * phi_no_bc);

    const auto [suppressed_full, suppressed_rhs] = MaxDeviationFromNoBCSolution(bc, true, assembly);
    EXPECT_LT(suppressed_full, 1.0e-6 * phi_no_bc);
    EXPECT_LT(suppressed_rhs, 1.0e-6 * phi_no_bc);
  }
}

} // namespace

TEST(DiffusionMIP, SuppressedBCsAddNoBoundaryTerms)
{
  CheckSuppressedBCs(Assembly::UNIT_CELL_MATRICES);
}

TEST(DiffusionMIP, SuppressedBCsAddNoBoundaryTermsWithQPoints)
{
  CheckSuppressedBCs(Assembly::QUADRATURE_POINTS);
}
