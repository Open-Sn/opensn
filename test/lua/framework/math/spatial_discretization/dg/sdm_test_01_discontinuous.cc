#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/lib/types.h"
#include "lua/lib/console.h"

using namespace opensn;

namespace unit_tests
{

/**Maps a face, in a discontinuous sense, using the spatial discretization.*/
int MapFaceNodeDisc(const CellMapping& cur_cell_mapping,
                    const CellMapping& adj_cell_mapping,
                    const std::vector<Vector3>& cc_node_locs,
                    const std::vector<Vector3>& ac_node_locs,
                    size_t ccf,
                    size_t acf,
                    size_t ccfi,
                    double epsilon = 1.0e-12);

double HPerpendicular(const CellMapping& cell_mapping, unsigned int f);

void
math_SDM_Test02_DisContinuous(const ParameterBlock& params)
{
  const double penalty_factor =
    params.Has("penalty_factor") ? params.GetParamValue<double>("penalty_factor") : 4.0;

  const bool export_vtk = params.Has("export_vtk") and params.GetParamValue<bool>("export_vtk");

  // Get grid
  auto grid = params.GetParamValue<std::shared_ptr<MeshContinuum>>("mesh");

  opensn::log.Log() << "Global num cells: " << grid->GetGlobalNumberOfCells();

  // Make SDM method
  const auto sdm_type = params.GetParamValue<std::string>("sdm_type");

  std::shared_ptr<SpatialDiscretization> sdm_ptr;

  {
    using namespace opensn;
    if (sdm_type == "PWLD")
      sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
    else
      OpenSnInvalidArgument("Unsupported sdm_type \"" + sdm_type + "\"");
  }

  auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_global_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num globl DOFs: " << num_global_dofs;

  // Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_global_dofs);
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
  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricFiniteElementData();
    const size_t num_nodes = cell_mapping.GetNumNodes();

    const auto& cc_nodes = cell_mapping.GetNodeLocations();

    const double D = 1.0;

    const auto [domain_nodes, bndry_nodes] = sdm.MakeCellInternalAndBndryNodeIDs(cell);

    DenseMatrix<double> Acell(num_nodes, num_nodes, 0.0);
    Vector<double> cell_rhs(num_nodes, 0.0);

    // Assemble continuous kernels
    {
      const auto& shape = qp_data.GetShapeValues();
      const auto& shape_grad = qp_data.GetShapeGradValues();
      const auto& JxW = qp_data.GetJxWValues();
      for (size_t i = 0; i < num_nodes; ++i)
      {
        if (bndry_nodes.find(i) != bndry_nodes.end())
          continue;
        for (size_t j = 0; j < num_nodes; ++j)
        {
          if (bndry_nodes.find(j) != bndry_nodes.end())
            continue;
          double entry_aij = 0.0;
          for (size_t qp : qp_data.GetQuadraturePointIndices())
            entry_aij += shape_grad[i][qp].Dot(shape_grad[j][qp]) * JxW[qp];

          Acell(i, j) = entry_aij;
        } // for j
        for (size_t qp : qp_data.GetQuadraturePointIndices())
          cell_rhs(i) += 1.0 * shape[i][qp] * JxW[qp];
      } // for i
    }   // continuous kernels

    const size_t num_faces = cell.faces.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const auto& n_f = face.normal;
      const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      const auto fqp_data = cell_mapping.MakeSurfaceFiniteElementData(f);

      const double hm = HPerpendicular(cell_mapping, f);

      if (face.has_neighbor)
      {
        const auto& adj_cell = grid->cells[face.neighbor_id];
        const auto& adj_cell_mapping = sdm.GetCellMapping(adj_cell);
        const auto ac_nodes = adj_cell_mapping.GetNodeLocations();
        const size_t acf = MeshContinuum::MapCellFace(cell, adj_cell, f);
        const double hp = HPerpendicular(adj_cell_mapping, acf);

        // Compute kappa
        double kappa = 1.0;
        if (cell.GetType() == CellType::SLAB)
          kappa = 2.0 * penalty_factor * (D / hp + D / hm) * 0.5;
        if (cell.GetType() == CellType::POLYGON)
          kappa = 2.0 * penalty_factor * (D / hp + D / hm) * 0.5;
        if (cell.GetType() == CellType::POLYHEDRON)
          kappa = 4.0 * penalty_factor * (D / hp + D / hm) * 0.5;

        // Assembly penalty terms
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
            for (size_t qp : fqp_data.GetQuadraturePointIndices())
              aij +=
                kappa * fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) * fqp_data.JxW(qp);

            MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A, imap, jpmap, -aij, ADD_VALUES);
          } // for fj
        }   // for fi

        // Assemble gradient terms
        // For the following comments we use the notation:
        // Dk = 0.5* n dot nabla bk

        // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
        for (int i = 0; i < num_nodes; ++i)
        {
          const int64_t imap = sdm.MapDOF(cell, i, OneDofPerNode, 0, 0);

          for (int fj = 0; fj < num_face_nodes; ++fj)
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

            Vector3 vec_aij;
            for (size_t qp : fqp_data.GetQuadraturePointIndices())
              vec_aij += fqp_data.ShapeValue(jm, qp) * fqp_data.ShapeGrad(i, qp) * fqp_data.JxW(qp);
            const double aij = -0.5 * D * n_f.Dot(vec_aij);

            MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A, imap, jpmap, -aij, ADD_VALUES);
          } // for fj
        }   // for i

        // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
        for (int fi = 0; fi < num_face_nodes; ++fi)
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

          for (int j = 0; j < num_nodes; ++j)
          {
            const int64_t jmap = sdm.MapDOF(cell, j, OneDofPerNode, 0, 0);

            Vector3 vec_aij;
            for (size_t qp : fqp_data.GetQuadraturePointIndices())
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

        // Compute kappa
        double kappa = 1.0;
        if (cell.GetType() == CellType::SLAB)
          kappa = 4.0 * penalty_factor * D / hm;
        if (cell.GetType() == CellType::POLYGON)
          kappa = 4.0 * penalty_factor * D / hm;
        if (cell.GetType() == CellType::POLYHEDRON)
          kappa = 8.0 * penalty_factor * D / hm;

        // Assembly penalty terms
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const int64_t imap = sdm.MapDOF(cell, i, OneDofPerNode, 0, 0);

          for (size_t fj = 0; fj < num_face_nodes; ++fj)
          {
            const int jm = cell_mapping.MapFaceNode(f, fj);
            const int64_t jmmap = sdm.MapDOF(cell, jm, OneDofPerNode, 0, 0);

            double aij = 0.0;
            for (size_t qp : fqp_data.GetQuadraturePointIndices())
              aij +=
                kappa * fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) * fqp_data.JxW(qp);
            double aij_bc_value = aij * bc_value;

            MatSetValue(A, imap, jmmap, aij, ADD_VALUES);
            VecSetValue(b, imap, aij_bc_value, ADD_VALUES);
          } // for fj
        }   // for fi

        // Assemble gradient terms
        // For the following comments we use the notation:
        // Dk = n dot nabla bk

        // D* n dot (b_j^+ - b_j^-)*nabla b_i^-
        for (size_t i = 0; i < num_nodes; ++i)
        {
          const int64_t imap = sdm.MapDOF(cell, i, OneDofPerNode, 0, 0);

          for (size_t j = 0; j < num_nodes; ++j)
          {
            const int64_t jmap = sdm.MapDOF(cell, j, OneDofPerNode, 0, 0);

            Vector3 vec_aij;
            for (size_t qp : fqp_data.GetQuadraturePointIndices())
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

    // Develop node mapping
    std::vector<int64_t> imap(num_nodes, 0); // node-mapping
    for (size_t i = 0; i < num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);

    // Assembly into system
    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t j = 0; j < num_nodes; ++j)
        MatSetValue(A, imap[i], imap[j], Acell(i, j), ADD_VALUES);

      VecSetValue(b, imap[i], cell_rhs(i), ADD_VALUES);
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
    CreateCommonKrylovSolverSetup(A, "PWLCDiffSolver", KSPCG, PCHYPRE, 0.0, 1.0e-9, 100);

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

  if (grid->GetDimension() == 2)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.6");
  if (grid->GetDimension() == 3)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.7");

  for (const auto& option : pc_options)
    PetscOptionsInsertString(nullptr, ("-" + option).c_str());

  PCSetFromOptions(pc);
  KSPSetFromOptions(petsc_solver.ksp);

  // Solve
  KSPSolve(petsc_solver.ksp, b, x);
  KSPConvergedReason reason;
  KSPGetConvergedReason(petsc_solver.ksp, &reason);
  if (reason == KSP_CONVERGED_ATOL)
    opensn::log.Log() << "Converged";

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
}

int
MapFaceNodeDisc(const CellMapping& cur_cell_mapping,
                const CellMapping& adj_cell_mapping,
                const std::vector<Vector3>& cc_node_locs,
                const std::vector<Vector3>& ac_node_locs,
                size_t ccf,
                size_t acf,
                size_t ccfi,
                double epsilon)
{
  const int i = cur_cell_mapping.MapFaceNode(ccf, ccfi);
  const auto& node_i_loc = cc_node_locs[i];

  const size_t adj_face_num_nodes = adj_cell_mapping.GetNumFaceNodes(acf);

  for (size_t fj = 0; fj < adj_face_num_nodes; ++fj)
  {
    const int j = adj_cell_mapping.MapFaceNode(acf, fj);
    if ((node_i_loc - ac_node_locs[j]).NormSquare() < epsilon)
      return j;
  }

  throw std::logic_error(
    "lbs::acceleration::DiffusionMIPSolver::MapFaceNodeDisc: Mapping failure.");
}

double
HPerpendicular(const CellMapping& cell_mapping, unsigned int f)
{
  const auto& cell = cell_mapping.GetCell();
  double hp;

  const auto num_faces = cell.faces.size();
  const auto num_vertices = cell.vertex_ids.size();

  const auto volume = cell.volume;
  const auto face_area = cell.faces.at(f).area;

  /**Lambda to compute surface area.*/
  auto ComputeSurfaceArea = [&cell, &num_faces]()
  {
    double surface_area = 0.0;
    for (size_t fr = 0; fr < num_faces; ++fr)
      surface_area += cell.faces[fr].area;

    return surface_area;
  };

  if (cell.GetType() == CellType::SLAB)
  {
    hp = volume / 2.0;
  }
  else if (cell.GetType() == CellType::POLYGON)
  {
    if (num_faces == 3)
      hp = 2.0 * volume / face_area;
    else if (num_faces == 4)
      hp = volume / face_area;
    else // Nv > 4
    {
      const auto surface_area = ComputeSurfaceArea();

      if (num_faces % 2 == 0)
        hp = 4.0 * volume / surface_area;
      else
      {
        hp = 2.0 * volume / surface_area;
        hp +=
          sqrt(2.0 * volume /
               (static_cast<double>(num_faces) * sin(2.0 * M_PI / static_cast<double>(num_faces))));
      }
    }
  }
  else if (cell.GetType() == CellType::POLYHEDRON)
  {
    const auto surface_area = ComputeSurfaceArea();

    if (num_faces == 4) // Tet
      hp = 3 * volume / surface_area;
    else if (num_faces == 6 and num_vertices == 8) // Hex
      hp = volume / surface_area;
    else // Polyhedron
      hp = 6 * volume / surface_area;
  } // Polyhedron
  else
    throw std::logic_error("lbs::acceleration::DiffusionMIPSolver::HPerpendicular: "
                           "Unsupported cell type in call to HPerpendicular");

  return hp;
}

BIND_FUNCTION(unit_tests, math_SDM_Test02_DisContinuous);

} //  namespace unit_tests
