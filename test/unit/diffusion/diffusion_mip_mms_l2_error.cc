// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "gmock/gmock.h"
#include "test/unit/common/mesh_builders.h"

#include "modules/diffusion/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"

#include "framework/field_functions/field_function_grid_based.h"
#include "framework/logging/log.h"
#include "framework/math/functions/function.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

double
MMS_phi(const Vector3& pt)
{
  return std::cos(M_PI * pt.x) + std::cos(M_PI * pt.y);
}

double
MMS_q(const Vector3& pt)
{
  return M_PI * M_PI * (std::cos(M_PI * pt.x) + std::cos(M_PI * pt.y));
}

std::pair<double, double>
SimTest_IP_MMS_L2error(std::shared_ptr<MeshContinuum> grid)
{
  using MatID2XSMap = std::map<unsigned int, Multigroup_D_and_sigR>;

  opensn::log.Log() << "SimTest_IP_MMS_L2error";
  opensn::log.Log() << "Global num cells: " << grid->GetGlobalNumberOfCells();

  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
  const auto& sdm = *sdm_ptr;
  const auto& one_dof_per_node = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(one_dof_per_node);
  const size_t num_global_dofs = sdm.GetNumGlobalDOFs(one_dof_per_node);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num global DOFs: " << num_global_dofs;

  std::map<uint64_t, BoundaryCondition> bcs;
  bcs[0] = {BCType::DIRICHLET, {2, 0, 0}}, bcs[1] = {BCType::DIRICHLET, {2, 0, 0}},
  bcs[2] = {BCType::DIRICHLET, {2, 0, 0}}, bcs[3] = {BCType::DIRICHLET, {2, 0, 0}},
  bcs[4] = {BCType::DIRICHLET, {2, 0, 0}}, bcs[5] = {BCType::DIRICHLET, {2, 0, 0}};

  MatID2XSMap matid_2_xs_map;
  matid_2_xs_map.insert(std::make_pair(0, Multigroup_D_and_sigR{{1.0}, {0.0}}));

  std::vector<UnitCellMatrices> unit_cell_matrices(grid->local_cells.size());
  for (const auto& cell : grid->local_cells)
    unit_cell_matrices[cell.local_id] = ComputeUnitCellIntegrals(sdm, cell);

  ScalarSpatialFunction mms_phi = MMS_phi;
  ScalarSpatialFunction mms_q = MMS_q;

  DiffusionMIPSolver solver("SimTest_IP_MMS_L2error",
                            sdm,
                            one_dof_per_node,
                            bcs,
                            matid_2_xs_map,
                            unit_cell_matrices,
                            false,
                            true);
  solver.options.verbose = true;
  solver.options.residual_tolerance = 1.0e-10;
  solver.options.perform_symmetry_check = true;
  solver.SetReferenceSolutionFunction(mms_phi);
  solver.SetSourceFunction(mms_q);
  solver.Initialize();

  std::vector<double> q_vector(num_local_dofs, 1.0);
  std::vector<double> x_vector(num_local_dofs, 0.0);
  solver.AssembleAand_b_wQpoints(q_vector);
  solver.Solve(x_vector);

  auto ff =
    std::make_shared<FieldFunctionGridBased>("Phi", sdm_ptr, one_dof_per_node.unknowns.front());

  auto compute_l2_error = [&grid, &sdm, &mms_phi, &ff](const std::vector<double>& solution)
  {
    ff->UpdateFieldVector(solution);
    const auto field_wg = ff->GetGhostedFieldVector();

    double local_error = 0.0;
    for (const auto& cell : grid->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.GetNumNodes();
      const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

      std::vector<double> nodal_phi(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; ++j)
      {
        const auto jmap = sdm.MapDOFLocal(cell, j);
        nodal_phi[j] = field_wg[jmap];
      }

      for (size_t qp : fe_vol_data.GetQuadraturePointIndices())
      {
        double phi_fem = 0.0;
        for (size_t j = 0; j < num_nodes; ++j)
          phi_fem += nodal_phi[j] * fe_vol_data.ShapeValue(j, qp);

        const double phi_true = mms_phi(fe_vol_data.QPointXYZ(qp));
        local_error += std::pow(phi_true - phi_fem, 2.0) * fe_vol_data.JxW(qp);
      }
    }

    double global_error = 0.0;
    opensn::mpi_comm.all_reduce(local_error, global_error, mpi::op::sum<double>());
    return std::sqrt(global_error);
  };

  const double full_assembly_error = compute_l2_error(x_vector);

  solver.Assemble_b_wQpoints(q_vector);
  solver.Solve(x_vector);
  const double rhs_assembly_error = compute_l2_error(x_vector);

  opensn::log.Log() << "Error: " << std::scientific << rhs_assembly_error
                    << " Num-cells: " << grid->GetGlobalNumberOfCells();

  return {full_assembly_error, rhs_assembly_error};
}

} // namespace

TEST(DiffusionMIP, MMSL2Convergence)
{
  if (opensn::mpi_comm.size() != 1)
    return;

  const std::array<unsigned int, 5> refinements = {5, 10, 20, 40, 80};
  const double length = 2.0;

  std::vector<double> errors;
  errors.reserve(refinements.size());

  for (const unsigned int num_cells : refinements)
  {
    auto grid = BuildSquareMesh(length, num_cells, -length / 2.0);
    grid->SetUniformBlockID(0);
    grid->SetOrthogonalBoundaries();

    const auto [full_assembly_error, rhs_assembly_error] = SimTest_IP_MMS_L2error(grid);
    errors.push_back(full_assembly_error);

    EXPECT_NEAR(rhs_assembly_error, full_assembly_error, 1.0e-12);
  }

  for (size_t k = 1; k < errors.size(); ++k)
  {
    EXPECT_LT(errors[k], errors[k - 1]);

    const double rate = std::log(errors[k - 1] / errors[k]) / std::log(2.0);
    opensn::log.Log() << "MIP MMS N=" << refinements[k] << " L2 rate=" << rate;
    EXPECT_GT(rate, 1.75);
  }

  opensn::mpi_comm.barrier();
}
