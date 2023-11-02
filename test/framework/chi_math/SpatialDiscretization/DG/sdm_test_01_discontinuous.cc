#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"
#include "math/SpatialDiscretization/FiniteElement/Lagrange/LagrangeDiscontinuous.h"
#include "math/PETScUtils/petsc_utils.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"

#define scdouble static_cast<double>

namespace chi_unit_tests
{

chi::InputParameters chi_math_SDM_Test02Syntax();
chi::ParameterBlock chi_math_SDM_Test02_DisContinuous(const chi::InputParameters& input_parameters);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chi_math_SDM_Test02_DisContinuous,
                        /*syntax_function=*/chi_math_SDM_Test02Syntax,
                        /*actual_function=*/chi_math_SDM_Test02_DisContinuous);

chi::InputParameters
chi_math_SDM_Test02Syntax()
{
  chi::InputParameters params;

  params.AddRequiredParameterBlock("arg0", "General parameters");

  return params;
}

int MapFaceNodeDisc(const chi_math::CellMapping& cur_cell_mapping,
                    const chi_math::CellMapping& adj_cell_mapping,
                    const std::vector<chi_mesh::Vector3>& cc_node_locs,
                    const std::vector<chi_mesh::Vector3>& ac_node_locs,
                    size_t ccf,
                    size_t acf,
                    size_t ccfi,
                    double epsilon = 1.0e-12);

double HPerpendicular(const chi_math::CellMapping& cell_mapping, unsigned int f);

chi::ParameterBlock
chi_math_SDM_Test02_DisContinuous(const chi::InputParameters& input_parameters)
{
  const chi::ParameterBlock& params = input_parameters.GetParam("arg0");

  const double penalty_factor =
    params.Has("penalty_factor") ? params.GetParamValue<double>("penalty_factor") : 4.0;

  const bool export_vtk = params.Has("export_vtk") && params.GetParamValue<bool>("export_vtk");

  //============================================= Get grid
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= Make SDM method
  const auto sdm_type = params.GetParamValue<std::string>("sdm_type");

  std::shared_ptr<chi_math::SpatialDiscretization> sdm_ptr;

  {
    using namespace chi_math::spatial_discretization;
    if (sdm_type == "PWLD") { sdm_ptr = PieceWiseLinearDiscontinuous::New(grid); }
    else if (sdm_type == "LagrangeD") { sdm_ptr = LagrangeDiscontinuous::New(grid); }
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
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  chi_math::PETScUtils::InitMatrixSparsity(A, nodal_nnz_in_diag, nodal_nnz_off_diag);

  //============================================= Assemble the system
  Chi::log.Log() << "Assembling system: ";
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto& cc_nodes = cell_mapping.GetNodeLocations();

    const double D = 1.0;

    const auto [domain_nodes, bndry_nodes] = sdm.MakeCellInternalAndBndryNodeIDs(cell);

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

    const size_t num_faces = cell.faces_.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      const auto& n_f = face.normal_;
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
      const auto fqp_data = cell_mapping.MakeSurfaceQuadraturePointData(f);

      const double hm = HPerpendicular(cell_mapping, f);

      typedef chi_mesh::MeshContinuum Grid;

      if (face.has_neighbor_)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id_];
        const auto& adj_cell_mapping = sdm.GetCellMapping(adj_cell);
        const auto ac_nodes = adj_cell_mapping.GetNodeLocations();
        const size_t acf = Grid::MapCellFace(cell, adj_cell, f);
        const double hp = HPerpendicular(adj_cell_mapping, acf);

        //========================= Compute kappa
        double kappa = 1.0;
        if (cell.Type() == chi_mesh::CellType::SLAB)
          kappa = 2.0 * penalty_factor * (D / hp + D / hm) * 0.5;
        if (cell.Type() == chi_mesh::CellType::POLYGON)
          kappa = 2.0 * penalty_factor * (D / hp + D / hm) * 0.5;
        if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
          kappa = 4.0 * penalty_factor * (D / hp + D / hm) * 0.5;

        //========================= Assembly penalty terms
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const int64_t imap = sdm.MapDOF(cell, i, OneDofPerNode, 0, 0);

          for (size_t fj = 0; fj < num_face_nodes; ++fj)
          {
            const int jm = cell_mapping.MapFaceNode(f, fj); // j-minus
            const int jp = MapFaceNodeDisc(cell_mapping,
                                           adj_cell_mapping,
                                           cc_nodes,
                                           ac_nodes,
                                           f,
                                           acf,
                                           fj); // j-plus
            const int64_t jmmap = sdm.MapDOF(cell, jm, OneDofPerNode, 0, 0);
            const int64_t jpmap = sdm.MapDOF(adj_cell, jp, OneDofPerNode, 0, 0);

            double aij = 0.0;
            for (size_t qp : fqp_data.QuadraturePointIndices())
              aij +=
                kappa * fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) * fqp_data.JxW(qp);

            MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A, imap, jpmap, -aij, ADD_VALUES);
          } // for fj
        }   // for fi

        //========================= Assemble gradient terms
        // For the following comments we use the notation:
        // Dk = 0.5* n dot nabla bk

        // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
        for (int i = 0; i < num_nodes; i++)
        {
          const int64_t imap = sdm.MapDOF(cell, i, OneDofPerNode, 0, 0);

          for (int fj = 0; fj < num_face_nodes; fj++)
          {
            const int jm = cell_mapping.MapFaceNode(f, fj); // j-minus
            const int jp = MapFaceNodeDisc(cell_mapping,
                                           adj_cell_mapping,
                                           cc_nodes,
                                           ac_nodes,
                                           f,
                                           acf,
                                           fj); // j-plus
            const int64_t jmmap = sdm.MapDOF(cell, jm, OneDofPerNode, 0, 0);
            const int64_t jpmap = sdm.MapDOF(adj_cell, jp, OneDofPerNode, 0, 0);

            chi_mesh::Vector3 vec_aij;
            for (size_t qp : fqp_data.QuadraturePointIndices())
              vec_aij += fqp_data.ShapeValue(jm, qp) * fqp_data.ShapeGrad(i, qp) * fqp_data.JxW(qp);
            const double aij = -0.5 * D * n_f.Dot(vec_aij);

            MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A, imap, jpmap, -aij, ADD_VALUES);
          } // for fj
        }   // for i

        // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
        for (int fi = 0; fi < num_face_nodes; fi++)
        {
          const int im = cell_mapping.MapFaceNode(f, fi); // i-minus
          const int ip = MapFaceNodeDisc(cell_mapping,
                                         adj_cell_mapping,
                                         cc_nodes,
                                         ac_nodes,
                                         f,
                                         acf,
                                         fi); // i-plus
          const int64_t immap = sdm.MapDOF(cell, im, OneDofPerNode, 0, 0);
          const int64_t ipmap = sdm.MapDOF(adj_cell, ip, OneDofPerNode, 0, 0);

          for (int j = 0; j < num_nodes; j++)
          {
            const int64_t jmap = sdm.MapDOF(cell, j, OneDofPerNode, 0, 0);

            chi_mesh::Vector3 vec_aij;
            for (size_t qp : fqp_data.QuadraturePointIndices())
              vec_aij += fqp_data.ShapeValue(im, qp) * fqp_data.ShapeGrad(j, qp) * fqp_data.JxW(qp);
            const double aij = -0.5 * D * n_f.Dot(vec_aij);

            MatSetValue(A, immap, jmap, aij, ADD_VALUES);
            MatSetValue(A, ipmap, jmap, -aij, ADD_VALUES);
          } // for j
        }   // for fi

      } // internal face
      else
      {
        const double bc_value = 0.0;

        //========================= Compute kappa
        double kappa = 1.0;
        if (cell.Type() == chi_mesh::CellType::SLAB) kappa = 4.0 * penalty_factor * D / hm;
        if (cell.Type() == chi_mesh::CellType::POLYGON) kappa = 4.0 * penalty_factor * D / hm;
        if (cell.Type() == chi_mesh::CellType::POLYHEDRON) kappa = 8.0 * penalty_factor * D / hm;

        //========================= Assembly penalty terms
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const int64_t imap = sdm.MapDOF(cell, i, OneDofPerNode, 0, 0);

          for (size_t fj = 0; fj < num_face_nodes; ++fj)
          {
            const int jm = cell_mapping.MapFaceNode(f, fj);
            const int64_t jmmap = sdm.MapDOF(cell, jm, OneDofPerNode, 0, 0);

            double aij = 0.0;
            for (size_t qp : fqp_data.QuadraturePointIndices())
              aij +=
                kappa * fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) * fqp_data.JxW(qp);
            double aij_bc_value = aij * bc_value;

            MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
            VecSetValue(b, imap, aij_bc_value, ADD_VALUES);
          } // for fj
        }   // for fi

        //========================= Assemble gradient terms
        // For the following comments we use the notation:
        // Dk = n dot nabla bk

        // D* n dot (b_j^+ - b_j^-)*nabla b_i^-
        for (size_t i = 0; i < num_nodes; i++)
        {
          const int64_t imap = sdm.MapDOF(cell, i, OneDofPerNode, 0, 0);

          for (size_t j = 0; j < num_nodes; j++)
          {
            const int64_t jmap = sdm.MapDOF(cell, j, OneDofPerNode, 0, 0);

            chi_mesh::Vector3 vec_aij;
            for (size_t qp : fqp_data.QuadraturePointIndices())
              vec_aij += fqp_data.ShapeValue(j, qp) * fqp_data.ShapeGrad(i, qp) * fqp_data.JxW(qp) +
                         fqp_data.ShapeValue(i, qp) * fqp_data.ShapeGrad(j, qp) * fqp_data.JxW(qp);
            const double aij = -D * n_f.Dot(vec_aij);

            double aij_bc_value = aij * bc_value;

            MatSetValue(A, imap, jmap, aij, ADD_VALUES);
            VecSetValue(b, imap, aij_bc_value, ADD_VALUES);
          } // for fj
        }   // for i
      }     // boundary face
    }       // for face

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
  auto petsc_solver =
    chi_math::PETScUtils::CreateCommonKrylovSolverSetup(A,                // Matrix
                                                        "PWLCDiffSolver", // Solver name
                                                        KSPCG,            // Solver type
                                                        PCHYPRE,          // Preconditioner type
                                                        1.0e-6, // Relative residual tolerance
                                                        1000);  // Max iterations

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

    chi_physics::FieldFunctionGridBased::ExportMultipleToVTK("ZSDM_Test", {ff});
  }

  return chi::ParameterBlock{};
}

// ###################################################################
/**Maps a face, in a discontinuous sense, using the spatial discretization.*/
int
MapFaceNodeDisc(const chi_math::CellMapping& cur_cell_mapping,
                const chi_math::CellMapping& adj_cell_mapping,
                const std::vector<chi_mesh::Vector3>& cc_node_locs,
                const std::vector<chi_mesh::Vector3>& ac_node_locs,
                size_t ccf,
                size_t acf,
                size_t ccfi,
                double epsilon /*=1.0e-12*/)
{
  const int i = cur_cell_mapping.MapFaceNode(ccf, ccfi);
  const auto& node_i_loc = cc_node_locs[i];

  const size_t adj_face_num_nodes = adj_cell_mapping.NumFaceNodes(acf);

  for (size_t fj = 0; fj < adj_face_num_nodes; ++fj)
  {
    const int j = adj_cell_mapping.MapFaceNode(acf, fj);
    if ((node_i_loc - ac_node_locs[j]).NormSquare() < epsilon) return j;
  }

  throw std::logic_error(
    "lbs::acceleration::DiffusionMIPSolver::MapFaceNodeDisc: Mapping failure.");
}

double
HPerpendicular(const chi_math::CellMapping& cell_mapping, unsigned int f)
{
  const auto& cell = cell_mapping.ReferenceCell();
  double hp;

  const size_t num_faces = cell.faces_.size();
  const size_t num_vertices = cell.vertex_ids_.size();

  const double volume = cell_mapping.CellVolume();
  const double face_area = cell_mapping.FaceArea(f);

  /**Lambda to compute surface area.*/
  auto ComputeSurfaceArea = [&cell_mapping, &num_faces]()
  {
    double surface_area = 0.0;
    for (size_t fr = 0; fr < num_faces; ++fr)
      surface_area += cell_mapping.FaceArea(fr);

    return surface_area;
  };

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell.Type() == chi_mesh::CellType::SLAB) hp = volume / 2.0;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    if (num_faces == 3) hp = 2.0 * volume / face_area;
    else if (num_faces == 4)
      hp = volume / face_area;
    else // Nv > 4
    {
      const double surface_area = ComputeSurfaceArea();

      if (num_faces % 2 == 0) hp = 4.0 * volume / surface_area;
      else
      {
        hp = 2.0 * volume / surface_area;
        hp += sqrt(2.0 * volume / (scdouble(num_faces) * sin(2.0 * M_PI / scdouble(num_faces))));
      }
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    const double surface_area = ComputeSurfaceArea();

    if (num_faces == 4) // Tet
      hp = 3 * volume / surface_area;
    else if (num_faces == 6 && num_vertices == 8) // Hex
      hp = volume / surface_area;
    else // Polyhedron
      hp = 6 * volume / surface_area;
  } // Polyhedron
  else
    throw std::logic_error("lbs::acceleration::DiffusionMIPSolver::HPerpendicular: "
                           "Unsupported cell type in call to HPerpendicular");

  return hp;
}

} // namespace chi_unit_tests