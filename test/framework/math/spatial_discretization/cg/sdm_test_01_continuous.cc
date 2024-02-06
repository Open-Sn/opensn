#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/math/spatial_discretization/finite_element/lagrange/lagrange_continuous.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace unit_tests
{

InputParameters math_SDM_Test01Syntax();
ParameterBlock math_SDM_Test01_Continuous(const InputParameters& input_parameters);

RegisterWrapperFunctionNamespace(unit_tests,
                                 math_SDM_Test01_Continuous,
                                 math_SDM_Test01Syntax,
                                 math_SDM_Test01_Continuous);

InputParameters
math_SDM_Test01Syntax()
{
  InputParameters params;

  params.AddRequiredParameterBlock("arg0", "General parameters");

  return params;
}

ParameterBlock
math_SDM_Test01_Continuous(const InputParameters& input_parameters)
{
  const ParameterBlock& params = input_parameters.GetParam("arg0");

  const bool export_vtk = params.Has("export_vtk") and params.GetParamValue<bool>("export_vtk");

  // Get grid
  auto grid_ptr = GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  opensn::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // Make SDM method
  const auto sdm_type = params.GetParamValue<std::string>("sdm_type");

  std::shared_ptr<SpatialDiscretization> sdm_ptr;
  bool is_DG = false;
  {
    using namespace opensn;
    if (sdm_type == "PWLC") sdm_ptr = PieceWiseLinearContinuous::New(grid);
    else if (sdm_type == "LagrangeC") { sdm_ptr = LagrangeContinuous::New(grid); }
    else
      ChiInvalidArgument("Unsupported sdm_type \"" + sdm_type + "\"");
  }

  auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  // Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_globl_dofs);
  Mat A;
  Vec x, b;

  A = CreateSquareMatrix(n, N);
  x = CreateVector(n, N);
  b = CreateVector(n, N);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  InitMatrixSparsity(A, nodal_nnz_in_diag, nodal_nnz_off_diag);

  // Assemble the system
  opensn::log.Log() << "Assembling system: ";
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto [domain_nodes, bndry_nodes] = sdm.MakeCellInternalAndBndryNodeIDs(cell);

    MatDbl Acell(num_nodes, VecDbl(num_nodes, 0.0));
    VecDbl cell_rhs(num_nodes, 0.0);

    // Assemble continuous kernels
    {
      const auto& shape = qp_data.ShapeValues();
      const auto& shape_grad = qp_data.ShapeGradValues();
      const auto& JxW = qp_data.JxW_Values();
      for (size_t i = 0; i < num_nodes; ++i)
      {
        if (bndry_nodes.find(i) != bndry_nodes.end()) continue;
        for (size_t j = 0; j < num_nodes; ++j)
        {
          if (bndry_nodes.find(j) != bndry_nodes.end()) continue;
          double entry_aij = 0.0;
          for (size_t qp : qp_data.QuadraturePointIndices())
            entry_aij += shape_grad[i][qp].Dot(shape_grad[j][qp]) * JxW[qp];

          Acell[i][j] = entry_aij;
        } // for j
        for (size_t qp : qp_data.QuadraturePointIndices())
          cell_rhs[i] += 1.0 * shape[i][qp] * JxW[qp];
      } // for i
    }   // continuous kernels

    // Apply dirichlet BCs
    for (auto i : bndry_nodes)
    {
      Acell[i][i] = 1.0;
      cell_rhs[i] = 0.0;
    }

    // Develop node mapping
    std::vector<int64_t> imap(num_nodes, 0); // node-mapping
    for (size_t i = 0; i < num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);

    // Assembly into system
    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t j = 0; j < num_nodes; ++j)
        MatSetValue(A, imap[i], imap[j], Acell[i][j], ADD_VALUES);

      VecSetValue(b, imap[i], cell_rhs[i], ADD_VALUES);
    } // for i
  }   // for cell

  opensn::log.Log() << "Global assembly";

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  opensn::log.Log() << "Done global assembly";

  // Create Krylov Solver
  opensn::log.Log() << "Solving: ";
  auto petsc_solver =
    CreateCommonKrylovSolverSetup(A, "PWLCDiffSolver", KSPCG, PCHYPRE, 1.0e-6, 1000);

  PC pc;
  KSPGetPC(petsc_solver.ksp, &pc);
  PCHYPRESetType(pc, "boomeramg");
  std::vector<std::string> pc_options = {"pc_hypre_boomeramg_agg_nl 1",
                                         "pc_hypre_boomeramg_P_max 4",
                                         "pc_hypre_boomeramg_grid_sweeps_coarse 1",
                                         "pc_hypre_boomeramg_max_levels 25",
                                         "pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi",
                                         "pc_hypre_boomeramg_coarsen_type HMIS",
                                         "pc_hypre_boomeramg_interp_type ext+i"};

  if (grid.Attributes() & DIMENSION_2)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.6");
  if (grid.Attributes() & DIMENSION_3)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.7");

  for (const auto& option : pc_options)
    PetscOptionsInsertString(nullptr, ("-" + option).c_str());

  PCSetFromOptions(pc);
  KSPSetFromOptions(petsc_solver.ksp);

  // Solve
  KSPSolve(petsc_solver.ksp, b, x);

  const char* reason;
  KSPGetConvergedReasonString(petsc_solver.ksp, &reason);
  opensn::log.Log() << "Done solving " << reason;

  // Extract PETSc vector
  std::vector<double> field;
  sdm.LocalizePETScVector(x, field, OneDofPerNode);

  double local_max = field.front();
  for (auto val : field)
    local_max = std::max(val, local_max);

  double global_max;
  opensn::mpi_comm.all_reduce(local_max, global_max, mpi::op::max<double>());

  opensn::log.Log() << "Nodal max = " << global_max;

  // Clean up
  KSPDestroy(&petsc_solver.ksp);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  opensn::log.Log() << "Done cleanup";

  // Create Field Function
  if (export_vtk)
  {
    auto ff =
      std::make_shared<FieldFunctionGridBased>("Phi", sdm_ptr, Unknown(UnknownType::SCALAR));

    ff->UpdateFieldVector(field);

    FieldFunctionGridBased::ExportMultipleToVTK("ZSDM_Test", {ff});
  }

  return ParameterBlock{};
}

} //  namespace unit_tests
