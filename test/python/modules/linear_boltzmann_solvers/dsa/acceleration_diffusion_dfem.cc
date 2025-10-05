// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/data_types/dense_matrix.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "test/python/src/bindings.h"

using namespace opensn;

namespace unit_tests
{

void
acceleration_Diffusion_DFEM(std::shared_ptr<MeshContinuum> grid)
{
  using MatID2XSMap = std::map<int, Multigroup_D_and_sigR>;
  opensn::log.Log() << "SimTest92_DSA";

  opensn::log.Log() << "Global num cells: " << grid->GetGlobalNumberOfCells();

  // Make SDM
  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_global_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num globl DOFs: " << num_global_dofs;

  // Make Boundary conditions
  std::map<uint64_t, BoundaryCondition> bcs;
  bcs[0] = {BCType::DIRICHLET, {2, 0, 0}}, bcs[1] = {BCType::DIRICHLET, {2, 0, 0}},
  bcs[2] = {BCType::DIRICHLET, {2, 0, 0}}, bcs[3] = {BCType::DIRICHLET, {2, 0, 0}},
  bcs[4] = {BCType::DIRICHLET, {2, 0, 0}}, bcs[5] = {BCType::DIRICHLET, {2, 0, 0}};

  MatID2XSMap matid_2_xs_map;
  matid_2_xs_map.insert(std::make_pair(0, Multigroup_D_and_sigR{{1.0}, {0.0}}));

  std::vector<UnitCellMatrices> unit_cell_matrices;
  unit_cell_matrices.resize(grid->local_cells.size());

  // Build unit integrals
  for (const auto& cell : grid->local_cells)
  {
    unit_cell_matrices[cell.local_id] = ComputeUnitCellIntegrals(sdm, cell);
  }

  // Retrieve the functions defined in the global namespace.
  py::module main_module = py::module::import("__main__");
  py::object mms_phi_fn = main_module.attr("mms_phi_fn");
  py::object mms_q_fn = main_module.attr("mms_q_fn");
  ScalarSpatialFunction mms_phi = mms_phi_fn.cast<ScalarSpatialFunction>();
  ScalarSpatialFunction mms_q = mms_q_fn.cast<ScalarSpatialFunction>();

  // Make solver
  DiffusionMIPSolver solver(
    "SimTest92_DSA", sdm, OneDofPerNode, bcs, matid_2_xs_map, unit_cell_matrices, false, true);
  solver.options.verbose = true;
  solver.options.residual_tolerance = 1.0e-10;
  solver.options.perform_symmetry_check = true;
  solver.SetReferenceSolutionFunction(mms_phi);
  solver.SetSourceFunction(mms_q);
  solver.Initialize();
  opensn::log.Log() << "Done constructing solver" << std::endl;

  // Assemble and solve
  std::vector<double> q_vector(num_local_dofs, 1.0);
  std::vector<double> x_vector(num_local_dofs, 0.0);
  solver.AssembleAand_b_wQpoints(q_vector);
  solver.Solve(x_vector);

  // Assemble and solver again
  solver.Assemble_b_wQpoints(q_vector);
  solver.Solve(x_vector);

  // Make Field-Function
  auto ff =
    std::make_shared<FieldFunctionGridBased>("Phi", sdm_ptr, OneDofPerNode.unknowns.front());
  ff->UpdateFieldVector(x_vector);
  FieldFunctionGridBased::ExportMultipleToPVTU("SimTest_92a_DSA", {ff});

  // Compute error
  // First get ghosted values
  const auto field_wg = ff->GetGhostedFieldVector();
  double local_error = 0.0;
  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    // Grab nodal phi values
    std::vector<double> nodal_phi(num_nodes, 0.0);
    for (size_t j = 0; j < num_nodes; ++j)
    {
      const auto jmap = sdm.MapDOFLocal(cell, j);
      nodal_phi[j] = field_wg[jmap];
    } // for j

    // Quadrature loop
    for (size_t qp : fe_vol_data.GetQuadraturePointIndices())
    {
      double phi_fem = 0.0;
      for (size_t j = 0; j < num_nodes; ++j)
        phi_fem += nodal_phi[j] * fe_vol_data.ShapeValue(j, qp);

      double phi_true = mms_phi(fe_vol_data.QPointXYZ(qp));

      local_error += std::pow(phi_true - phi_fem, 2.0) * fe_vol_data.JxW(qp);
    }
  } // for cell

  double global_error = 0.0;
  opensn::mpi_comm.all_reduce(local_error, global_error, mpi::op::sum<double>());

  global_error = std::sqrt(global_error);

  opensn::log.Log() << "Error: " << std::scientific << global_error
                    << " Num-cells: " << grid->GetGlobalNumberOfCells();
}

} //  namespace unit_tests
