#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"
#include "math/PETScUtils/petsc_utils.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"

#include "chi_lua.h"

namespace chi_unit_sim_tests
{

chi::ParameterBlock
chiSimTest04_PWLC(const chi::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chiSimTest04_PWLC,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/chiSimTest04_PWLC);

/**This is a simple test of the Finite Volume spatial discretization applied
 * to Laplace's problem but with a manufactured solution. */
chi::ParameterBlock
chiSimTest04_PWLC(const chi::InputParameters& params)
{
  Chi::log.Log() << "Coding Tutorial 4";

  //============================================= Get grid
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= Make SDM
  typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;
  SDMPtr sdm_ptr = chi_math::spatial_discretization::PieceWiseLinearContinuous::New(grid);
  const auto& sdm = *sdm_ptr;

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

  //============================================= Source lambda
  lua_State* L = chi::Console::GetInstance().GetConsoleState();
  auto CallLuaXYZFunction =
    [&L](const std::string& lua_func_name, const chi_mesh::Vector3& xyz)
  {
    //============= Load lua function
    lua_getglobal(L, lua_func_name.c_str());

    //============= Error check lua function
    if (not lua_isfunction(L, -1))
      throw std::logic_error(
        "CallLuaXYZFunction attempted to access lua-function, " +
        lua_func_name +
        ", but it seems the function"
        " could not be retrieved.");

    //============= Push arguments
    lua_pushnumber(L, xyz.x);
    lua_pushnumber(L, xyz.y);
    lua_pushnumber(L, xyz.z);

    //============= Call lua function
    // 3 arguments, 1 result (double), 0=original error object
    double lua_return = 0.0;
    if (lua_pcall(L, 3, 1, 0) == 0)
    {
      LuaCheckNumberValue("CallLuaXYZFunction", L, -1);
      lua_return = lua_tonumber(L, -1);
    }
    else
      throw std::logic_error(
        "CallLuaXYZFunction attempted to call lua-function, " + lua_func_name +
        ", but the call failed." + xyz.PrintStr());

    lua_pop(L, 1); // pop the double, or error code

    return lua_return;
  };

  //============================================= Assemble the system
  Chi::log.Log() << "Assembling system: ";
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const auto cell_node_xyzs = cell_mapping.GetNodeLocations();

    const size_t num_nodes = cell_mapping.NumNodes();
    MatDbl Acell(num_nodes, VecDbl(num_nodes, 0.0));
    VecDbl cell_rhs(num_nodes, 0.0);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t j = 0; j < num_nodes; ++j)
      {
        double entry_aij = 0.0;
        for (size_t qp : qp_data.QuadraturePointIndices())
        {
          entry_aij += qp_data.ShapeGrad(i, qp).Dot(qp_data.ShapeGrad(j, qp)) *
                       qp_data.JxW(qp);
        } // for qp
        Acell[i][j] = entry_aij;
      } // for j
      for (size_t qp : qp_data.QuadraturePointIndices())
        cell_rhs[i] += CallLuaXYZFunction("MMS_q", qp_data.QPointXYZ(qp)) *
                       qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
    } // for i

    //======================= Flag nodes for being on dirichlet boundary
    std::vector<bool> node_boundary_flag(num_nodes, false);
    const size_t num_faces = cell.faces_.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      if (face.has_neighbor_) continue;

      const size_t num_face_nodes = face.vertex_ids_.size();
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const uint i = cell_mapping.MapFaceNode(f, fi);
        node_boundary_flag[i] = true;
      } // for fi
    }   // for face f

    //======================= Develop node mapping
    std::vector<int64_t> imap(num_nodes, 0); // node-mapping
    for (size_t i = 0; i < num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);

    //======================= Assembly into system
    for (size_t i = 0; i < num_nodes; ++i)
    {
      if (node_boundary_flag[i]) // if dirichlet node
      {
        MatSetValue(A, imap[i], imap[i], 1.0, ADD_VALUES);
        double bval = CallLuaXYZFunction("MMS_phi", cell_node_xyzs[i]);
        VecSetValue(b, imap[i], bval, ADD_VALUES);
      }
      else
      {
        for (size_t j = 0; j < num_nodes; ++j)
        {
          if (not node_boundary_flag[j])
            MatSetValue(A, imap[i], imap[j], Acell[i][j], ADD_VALUES);
          else
          {
            double bval = CallLuaXYZFunction("MMS_phi", cell_node_xyzs[j]);
            VecSetValue(b, imap[i], -Acell[i][j] * bval, ADD_VALUES);
          }
        } // for j
        VecSetValue(b, imap[i], cell_rhs[i], ADD_VALUES);
      }
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
    PCGAMG,           // Preconditioner type
    1.0e-6,           // Relative residual tolerance
    1000);            // Max iterations

  //============================================= Solve
  KSPSolve(petsc_solver.ksp, b, x);

  Chi::log.Log() << "Done solving";

  //============================================= Extract PETSc vector
  std::vector<double> field;
  sdm.LocalizePETScVector(x, field, OneDofPerNode);

  //============================================= Clean up
  KSPDestroy(&petsc_solver.ksp);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  Chi::log.Log() << "Done cleanup";

  //============================================= Create Field Function
  auto ff = std::make_shared<chi_physics::FieldFunctionGridBased>(
    "Phi",                                           // Text name
    sdm_ptr,                                         // Spatial Discr.
    chi_math::Unknown(chi_math::UnknownType::SCALAR) // Unknown
  );

  ff->UpdateFieldVector(field);

  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK("CodeTut4_PWLC",
                                                           {ff});

  //============================================= Compute error
  // First get ghosted values
  const auto field_wg = ff->GetGhostedFieldVector();

  double local_error = 0.0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();

    //======================= Grab nodal phi values
    std::vector<double> nodal_phi(num_nodes, 0.0);
    for (size_t j = 0; j < num_nodes; ++j)
    {
      const int64_t jmap = sdm.MapDOFLocal(cell, j);
      nodal_phi[j] = field_wg[jmap];
    } // for j

    //======================= Quadrature loop
    for (size_t qp : qp_data.QuadraturePointIndices())
    {
      double phi_fem = 0.0;
      for (size_t j = 0; j < num_nodes; ++j)
        phi_fem += nodal_phi[j] * qp_data.ShapeValue(j, qp);

      double phi_true = CallLuaXYZFunction("MMS_phi", qp_data.QPointXYZ(qp));

      local_error += std::pow(phi_true - phi_fem, 2.0) * qp_data.JxW(qp);
    }
  } // for cell

  double global_error = 0.0;
  MPI_Allreduce(&local_error,  // sendbuf
                &global_error, // recvbuf
                1,
                MPI_DOUBLE,      // count+datatype
                MPI_SUM,         // operation
                Chi::mpi.comm); // communicator

  global_error = std::sqrt(global_error);

  Chi::log.Log() << "Error: " << std::scientific << global_error
                 << " Num-cells: " << grid.GetGlobalNumberOfCells();

  return chi::ParameterBlock();
}

} // namespace chi_unit_sim_tests
