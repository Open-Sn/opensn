#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/acceleration/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_structs.h"

#include "framework/physics/field_function/field_function_grid_based.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "lua/framework/console/console.h"
#include "lua/framework/math/functions/lua_scalar_spatial_function.h"

using namespace opensn;
using namespace opensnlua;

namespace unit_sim_tests
{

namespace
{

std::shared_ptr<LuaScalarSpatialFunction>
CreateFunction(const std::string& function_name)
{
  ParameterBlock blk;
  blk.AddParameter("lua_function_name", function_name);
  InputParameters params = LuaScalarSpatialFunction::GetInputParameters();
  params.AssignParameters(blk);
  return std::make_shared<LuaScalarSpatialFunction>(params);
}

} // namespace

ParameterBlock acceleration_Diffusion_DFEM(const InputParameters& params);

RegisterWrapperFunction(chi_unit_tests,
                        acceleration_Diffusion_DFEM,
                        nullptr,
                        acceleration_Diffusion_DFEM);

ParameterBlock
acceleration_Diffusion_DFEM(const InputParameters&)
{
  typedef std::map<int, lbs::Multigroup_D_and_sigR> MatID2XSMap;
  opensn::log.Log() << "chiSimTest92_DSA";

  // Get grid
  auto grid_ptr = GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  opensn::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // Make SDM
  typedef std::shared_ptr<SpatialDiscretization> SDMPtr;
  SDMPtr sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  // Make Boundary conditions
  typedef lbs::BoundaryCondition BC;
  std::map<uint64_t, BC> bcs;
  bcs[0] = {lbs::BCType::DIRICHLET, {2, 0, 0}}, bcs[1] = {lbs::BCType::DIRICHLET, {2, 0, 0}},
  bcs[2] = {lbs::BCType::DIRICHLET, {2, 0, 0}}, bcs[3] = {lbs::BCType::DIRICHLET, {2, 0, 0}},
  bcs[4] = {lbs::BCType::DIRICHLET, {2, 0, 0}}, bcs[5] = {lbs::BCType::DIRICHLET, {2, 0, 0}};

  MatID2XSMap matid_2_xs_map;
  matid_2_xs_map.insert(std::make_pair(0, lbs::Multigroup_D_and_sigR{{1.0}, {0.0}}));

  std::vector<lbs::UnitCellMatrices> unit_cell_matrices;
  unit_cell_matrices.resize(grid.local_cells.size());

  // Build unit integrals
  typedef std::vector<Vector3> VecVec3;
  typedef std::vector<VecVec3> MatVec3;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t cell_num_faces = cell.faces_.size();
    const size_t cell_num_nodes = cell_mapping.NumNodes();
    const auto vol_qp_data = cell_mapping.MakeVolumetricQuadraturePointData();

    MatDbl IntV_gradshapeI_gradshapeJ(cell_num_nodes, VecDbl(cell_num_nodes));
    MatDbl IntV_shapeI_shapeJ(cell_num_nodes, VecDbl(cell_num_nodes));
    VecDbl IntV_shapeI(cell_num_nodes);

    std::vector<MatDbl> IntS_shapeI_shapeJ(cell_num_faces);
    std::vector<MatVec3> IntS_shapeI_gradshapeJ(cell_num_faces);
    std::vector<VecDbl> IntS_shapeI(cell_num_faces);

    // Volume integrals
    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : vol_qp_data.QuadraturePointIndices())
        {
          IntV_gradshapeI_gradshapeJ[i][j] +=
            vol_qp_data.ShapeGrad(i, qp).Dot(vol_qp_data.ShapeGrad(j, qp)) *
            vol_qp_data.JxW(qp); // K-matrix

          IntV_shapeI_shapeJ[i][j] += vol_qp_data.ShapeValue(i, qp) *
                                      vol_qp_data.ShapeValue(j, qp) *
                                      vol_qp_data.JxW(qp); // M-matrix
        }                                                  // for qp
      }                                                    // for j

      for (const auto& qp : vol_qp_data.QuadraturePointIndices())
      {
        IntV_shapeI[i] += vol_qp_data.ShapeValue(i, qp) * vol_qp_data.JxW(qp);
      } // for qp
    }   // for i

    //  surface integrals
    for (size_t f = 0; f < cell_num_faces; ++f)
    {
      const auto faces_qp_data = cell_mapping.MakeSurfaceQuadraturePointData(f);
      IntS_shapeI_shapeJ[f].resize(cell_num_nodes, VecDbl(cell_num_nodes));
      IntS_shapeI[f].resize(cell_num_nodes);
      IntS_shapeI_gradshapeJ[f].resize(cell_num_nodes, VecVec3(cell_num_nodes));

      for (unsigned int i = 0; i < cell_num_nodes; ++i)
      {
        for (unsigned int j = 0; j < cell_num_nodes; ++j)
        {
          for (const auto& qp : faces_qp_data.QuadraturePointIndices())
          {
            IntS_shapeI_shapeJ[f][i][j] += faces_qp_data.ShapeValue(i, qp) *
                                           faces_qp_data.ShapeValue(j, qp) * faces_qp_data.JxW(qp);
            IntS_shapeI_gradshapeJ[f][i][j] += faces_qp_data.ShapeValue(i, qp) *
                                               faces_qp_data.ShapeGrad(j, qp) *
                                               faces_qp_data.JxW(qp);
          } // for qp
        }   // for j

        for (const auto& qp : faces_qp_data.QuadraturePointIndices())
        {
          IntS_shapeI[f][i] += faces_qp_data.ShapeValue(i, qp) * faces_qp_data.JxW(qp);
        } // for qp
      }   // for i
    }     // for f

    unit_cell_matrices[cell.local_id_] = lbs::UnitCellMatrices{IntV_gradshapeI_gradshapeJ,
                                                               {},
                                                               IntV_shapeI_shapeJ,
                                                               IntV_shapeI,

                                                               IntS_shapeI_shapeJ,
                                                               IntS_shapeI_gradshapeJ,
                                                               IntS_shapeI};
  } // for cell

  auto mms_phi_function = CreateFunction("MMS_phi");
  opensn::function_stack.push_back(mms_phi_function);

  auto mms_q_function = CreateFunction("MMS_q");
  opensn::function_stack.push_back(mms_q_function);

  // Make solver
  lbs::DiffusionMIPSolver solver(
    "SimTest92_DSA", sdm, OneDofPerNode, bcs, matid_2_xs_map, unit_cell_matrices, true);
  solver.options.verbose = true;
  solver.options.residual_tolerance = 1.0e-10;
  solver.options.perform_symmetry_check = true;

  solver.SetReferenceSolutionFunction(mms_phi_function);
  solver.SetSourceFunction(mms_q_function);

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
    std::make_shared<FieldFunctionGridBased>("Phi", sdm_ptr, OneDofPerNode.unknowns_.front());

  ff->UpdateFieldVector(x_vector);

  FieldFunctionGridBased::ExportMultipleToVTK("SimTest_92a_DSA", {ff});

  // Compute error
  // First get ghosted values
  const auto field_wg = ff->GetGhostedFieldVector();

  double local_error = 0.0;
  lua_State* L = Console::GetInstance().GetConsoleState();
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();

    // Grab nodal phi values
    std::vector<double> nodal_phi(num_nodes, 0.0);
    for (size_t j = 0; j < num_nodes; ++j)
    {
      const int64_t jmap = sdm.MapDOFLocal(cell, j);
      nodal_phi[j] = field_wg[jmap];
    } // for j

    // Quadrature loop
    for (size_t qp : qp_data.QuadraturePointIndices())
    {
      double phi_fem = 0.0;
      for (size_t j = 0; j < num_nodes; ++j)
        phi_fem += nodal_phi[j] * qp_data.ShapeValue(j, qp);

      double phi_true = mms_phi_function->Evaluate(qp_data.QPointXYZ(qp));

      local_error += std::pow(phi_true - phi_fem, 2.0) * qp_data.JxW(qp);
    }
  } // for cell

  double global_error = 0.0;
  MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, opensn::mpi.comm);

  global_error = std::sqrt(global_error);

  opensn::log.Log() << "Error: " << std::scientific << global_error
                    << " Num-cells: " << grid.GetGlobalNumberOfCells();

  return ParameterBlock();
}

} //  namespace unit_sim_tests
