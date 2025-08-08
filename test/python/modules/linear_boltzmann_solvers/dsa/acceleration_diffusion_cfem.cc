// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "modules/diffusion/diffusion_pwlc_solver.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "test/python/src/bindings.h"

using namespace opensn;

namespace unit_tests
{

void
acceleration_Diffusion_CFEM(std::shared_ptr<MeshContinuum> grid)
{
  using MatID2XSMap = std::map<int, Multigroup_D_and_sigR>;
  opensn::log.Log() << "SimTest92_DSA";

  opensn::log.Log() << "Global num cells: " << grid->GetGlobalNumberOfCells();

  // Make SDM
  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearContinuous::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalAndGhostDOFs(OneDofPerNode);
  const size_t num_global_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num globl DOFs: " << num_global_dofs;

  // Make Boundary conditions
  std::map<uint64_t, BoundaryCondition> bcs;
  bcs[0] = {BCType::ROBIN, {0.25, 0.5, 0}}, bcs[1] = {BCType::ROBIN, {0.25, 0.5, 0}},
  bcs[2] = {BCType::ROBIN, {0.25, 0.5, 0}}, bcs[3] = {BCType::ROBIN, {0.25, 0.5, 0}},
  bcs[4] = {BCType::ROBIN, {0.25, 0.5, 0}}, bcs[5] = {BCType::ROBIN, {0.25, 0.5, 0}};

  MatID2XSMap matid_2_xs_map;
  matid_2_xs_map.insert(std::make_pair(0, Multigroup_D_and_sigR{{1.0}, {0.0}}));

  std::vector<UnitCellMatrices> unit_cell_matrices;
  unit_cell_matrices.resize(grid->local_cells.size());

  // Build unit integrals
  for (const auto& cell : grid->local_cells)
  {
    unit_cell_matrices[cell.local_id] = ComputeUnitCellIntegrals(sdm, cell);
  }

  // Make solver
  DiffusionPWLCSolver solver("SimTest92b_DSA_PWLC",
                             sdm,
                             OneDofPerNode,
                             bcs,
                             matid_2_xs_map,
                             unit_cell_matrices,
                             false,
                             true);
  // TODO: For this to work, add MMS support into `lbs/acceleration/DiffusionSolver`
  // solver.options.ref_solution_function = "MMS_phi";
  // solver.options.source_function = "MMS_q";
  solver.options.verbose = true;
  solver.options.residual_tolerance = 1.0e-12;
  solver.options.perform_symmetry_check = true;

  solver.Initialize();

  opensn::log.Log() << "Done constructing solver" << std::endl;

  // Assemble and solve
  std::vector<double> q_vector(num_local_dofs, 1.0);
  std::vector<double> x_vector(num_local_dofs, 0.0);

  solver.AssembleAand_b(q_vector);
  solver.Solve(x_vector);

  // Assemble and solve again
  solver.Assemble_b(q_vector);
  solver.Solve(x_vector);

  // Make Field-Function
  auto ff =
    std::make_shared<FieldFunctionGridBased>("Phi", sdm_ptr, OneDofPerNode.unknowns.front());

  ff->UpdateFieldVector(x_vector);

  FieldFunctionGridBased::ExportMultipleToPVTU("SimTest_92b_DSA_PWLC", {ff});
}

} // namespace unit_tests
