#include "framework/mesh/MeshHandler/chi_meshhandler.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "framework/math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/acceleration.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/diffusion_pwlc_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_structs.h"

#include "framework/physics/FieldFunction/fieldfunction_gridbased.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "framework/console/chi_console.h"

namespace chi_unit_sim_tests
{

chi::ParameterBlock acceleration_Diffusion_CFEM(const chi::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/acceleration_Diffusion_CFEM,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/acceleration_Diffusion_CFEM);

chi::ParameterBlock
acceleration_Diffusion_CFEM(const chi::InputParameters&)
{
  typedef std::map<int, lbs::acceleration::Multigroup_D_and_sigR> MatID2XSMap;
  Chi::log.Log() << "chiSimTest92_DSA";

  //============================================= Get grid
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= Make SDM
  typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;
  SDMPtr sdm_ptr = chi_math::spatial_discretization::PieceWiseLinearContinuous::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalAndGhostDOFs(OneDofPerNode);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  Chi::log.Log() << "Num local DOFs: " << num_local_dofs;
  Chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  //============================================= Make Boundary conditions
  typedef lbs::acceleration::BoundaryCondition BC;
  std::map<uint64_t, BC> bcs;
  // bcs[0] = {lbs::acceleration::BCType::DIRICHLET,{1.0,0,0}},
  // bcs[1] = {lbs::acceleration::BCType::DIRICHLET,{1.0,0,0}},
  // bcs[2] = {lbs::acceleration::BCType::DIRICHLET,{1.0,0,0}},
  // bcs[3] = {lbs::acceleration::BCType::DIRICHLET,{1.0,0,0}},
  // bcs[4] = {lbs::acceleration::BCType::DIRICHLET,{1.0,0,0}},
  // bcs[5] = {lbs::acceleration::BCType::DIRICHLET,{1.0,0,0}};
  bcs[0] = {lbs::acceleration::BCType::ROBIN, {0.25, 0.5, 0}},
  bcs[1] = {lbs::acceleration::BCType::ROBIN, {0.25, 0.5, 0}},
  bcs[2] = {lbs::acceleration::BCType::ROBIN, {0.25, 0.5, 0}},
  bcs[3] = {lbs::acceleration::BCType::ROBIN, {0.25, 0.5, 0}},
  bcs[4] = {lbs::acceleration::BCType::ROBIN, {0.25, 0.5, 0}},
  bcs[5] = {lbs::acceleration::BCType::ROBIN, {0.25, 0.5, 0}};

  MatID2XSMap matid_2_xs_map;
  matid_2_xs_map.insert(std::make_pair(0, lbs::acceleration::Multigroup_D_and_sigR{{1.0}, {0.0}}));

  std::vector<lbs::UnitCellMatrices> unit_cell_matrices;
  unit_cell_matrices.resize(grid.local_cells.size());

  //============================================= Build unit integrals
  typedef std::vector<chi_mesh::Vector3> VecVec3;
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

    unit_cell_matrices[cell.local_id_] =
      lbs::UnitCellMatrices{IntV_gradshapeI_gradshapeJ, // K-matrix
                            {},                         // G-matrix
                            IntV_shapeI_shapeJ,         // M-matrix
                            IntV_shapeI,                // Vi-vectors

                            IntS_shapeI_shapeJ,     // face M-matrices
                            IntS_shapeI_gradshapeJ, // face G-matrices
                            IntS_shapeI};           // face Si-vectors
  }                                                 // for cell

  //============================================= Make solver
  lbs::acceleration::DiffusionPWLCSolver solver(
    "SimTest92b_DSA_PWLC", sdm, OneDofPerNode, bcs, matid_2_xs_map, unit_cell_matrices, true);
  solver.options.ref_solution_lua_function = "MMS_phi";
  solver.options.source_lua_function = "MMS_q";
  solver.options.verbose = true;
  solver.options.residual_tolerance = 1.0e-10;
  solver.options.perform_symmetry_check = true;

  solver.Initialize();

  Chi::log.Log() << "Done constructing solver" << std::endl;

  //============================================= Assemble and solve
  std::vector<double> q_vector(num_local_dofs, 1.0);
  std::vector<double> x_vector(num_local_dofs, 0.0);

  solver.AssembleAand_b(q_vector);
  solver.Solve(x_vector);

  //============================================= Assemble and solve again
  solver.Assemble_b(q_vector);
  solver.Solve(x_vector);

  //============================================= Make Field-Function
  auto ff = std::make_shared<chi_physics::FieldFunctionGridBased>(
    "Phi", sdm_ptr, OneDofPerNode.unknowns_.front());

  ff->UpdateFieldVector(x_vector);

  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK("SimTest_92b_DSA_PWLC", {ff});

  return chi::ParameterBlock();
}

} // namespace chi_unit_sim_tests
