#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"
#include "math/SpatialDiscretization/FiniteElement/Lagrange/LagrangeContinuous.h"
#include "math/PETScUtils/petsc_utils.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"

#define scdouble static_cast<double>

namespace chi_unit_tests
{

chi::InputParameters chi_math_SDM_Test01Syntax();
chi::ParameterBlock
chi_math_SDM_Test01_Continuous(const chi::InputParameters& input_parameters);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chi_math_SDM_Test01_Continuous,
                        /*syntax_function=*/chi_math_SDM_Test01Syntax,
                        /*actual_function=*/chi_math_SDM_Test01_Continuous);

chi::InputParameters chi_math_SDM_Test01Syntax()
{
  chi::InputParameters params;

  params.AddRequiredParameterBlock("arg0", "General parameters");

  return params;
}

chi::ParameterBlock
chi_math_SDM_Test01_Continuous(const chi::InputParameters& input_parameters)
{
  const chi::ParameterBlock& params = input_parameters.GetParam("arg0");

  const bool export_vtk =
    params.Has("export_vtk") && params.GetParamValue<bool>("export_vtk");

  //============================================= Get grid
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= Make SDM method
  const auto sdm_type = params.GetParamValue<std::string>("sdm_type");

  std::shared_ptr<chi_math::SpatialDiscretization> sdm_ptr;
  bool is_DG = false;
  {
    using namespace chi_math::spatial_discretization;
    if (sdm_type == "PWLC") sdm_ptr = PieceWiseLinearContinuous::New(grid);
    else if (sdm_type == "LagrangeC")
    {
      sdm_ptr = LagrangeContinuous::New(grid);
    }
    else
      ChiInvalidArgument("Unsupported sdm_type \"" + sdm_type + "\"");
  }

  auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  Chi::log.Log() << "Num local DOFs: " << num_local_dofs;
  Chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  //============================================= Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_globl_dofs);
  Mat A;
  Vec x, b;

  A = chi_math::PETScUtils::CreateSquareMatrix(n, N);
  x = chi_math::PETScUtils::CreateVector(n, N);
  b = chi_math::PETScUtils::CreateVector(n, N);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(
    nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  chi_math::PETScUtils::InitMatrixSparsity(
    A, nodal_nnz_in_diag, nodal_nnz_off_diag);

  //============================================= Assemble the system
  Chi::log.Log() << "Assembling system: ";
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto [domain_nodes, bndry_nodes] =
      sdm.MakeCellInternalAndBndryNodeIDs(cell);

    MatDbl Acell(num_nodes, VecDbl(num_nodes, 0.0));
    VecDbl cell_rhs(num_nodes, 0.0);

    //======================= Assemble continuous kernels
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

    //======================= Develop node mapping
    std::vector<int64_t> imap(num_nodes, 0); // node-mapping
    for (size_t i = 0; i < num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);

    //======================= Assembly into system
    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t j = 0; j < num_nodes; ++j)
        MatSetValue(A, imap[i], imap[j], Acell[i][j], ADD_VALUES);

      VecSetValue(b, imap[i], cell_rhs[i], ADD_VALUES);
    } // for i
  }   // for cell

  Chi::log.Log() << "Global assembly";

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  Chi::log.Log() << "Done global assembly";

  //============================================= Create Krylov Solver
  Chi::log.Log() << "Solving: ";
  auto petsc_solver = chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
    A,                // Matrix
    "PWLCDiffSolver", // Solver name
    KSPCG,            // Solver type
    PCHYPRE,           // Preconditioner type
    1.0e-6,           // Relative residual tolerance
    1000);            // Max iterations

  PC pc;
  KSPGetPC(petsc_solver.ksp, &pc);
  PCHYPRESetType(pc, "boomeramg");
  std::vector<std::string> pc_options = {
    "pc_hypre_boomeramg_agg_nl 1",
    "pc_hypre_boomeramg_P_max 4",
    "pc_hypre_boomeramg_grid_sweeps_coarse 1",
    "pc_hypre_boomeramg_max_levels 25",
    "pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi",
    "pc_hypre_boomeramg_coarsen_type HMIS",
    "pc_hypre_boomeramg_interp_type ext+i"};

  if (grid.Attributes() & chi_mesh::DIMENSION_2)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.6");
  if (grid.Attributes() & chi_mesh::DIMENSION_3)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.7");

  for (const auto& option : pc_options)
    PetscOptionsInsertString(nullptr, ("-" + option).c_str());

  PCSetFromOptions(pc);
  KSPSetFromOptions(petsc_solver.ksp);

  //============================================= Solve
  KSPSolve(petsc_solver.ksp, b, x);

  const char* reason;
  KSPGetConvergedReasonString(petsc_solver.ksp, &reason);
  Chi::log.Log() << "Done solving " << reason;


  //============================================= Extract PETSc vector
  std::vector<double> field;
  sdm.LocalizePETScVector(x, field, OneDofPerNode);

  double local_max = field.front();
  for (auto val : field)
    local_max = std::max(val, local_max);

  double global_max;
  MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, Chi::mpi.comm);

  Chi::log.Log() << "Nodal max = " << global_max;

  //============================================= Clean up
  KSPDestroy(&petsc_solver.ksp);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  Chi::log.Log() << "Done cleanup";

  //============================================= Create Field Function
  if (export_vtk)
  {
    auto ff = std::make_shared<chi_physics::FieldFunctionGridBased>(
      "Phi",                                           // Text name
      sdm_ptr,                                         // Spatial Discr.
      chi_math::Unknown(chi_math::UnknownType::SCALAR) // Unknown
    );

    ff->UpdateFieldVector(field);

    chi_physics::FieldFunctionGridBased::ExportMultipleToVTK(
      "ZSDM_Test", {ff});
  }

  return chi::ParameterBlock{};
}

} // namespace chi_unit_tests