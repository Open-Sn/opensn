#include "dfem_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "dfem_diffusion_bndry.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"

#define DefaultBCDirichlet BoundaryCondition{BCType::DIRICHLET,{0,0,0}}

//============================================= constructor
dfem_diffusion::Solver::Solver(const std::string& in_solver_name):
  chi_physics::Solver(in_solver_name, { {"max_iters", int64_t(500)   },
                                        {"residual_tolerance", 1.0e-2}})
{}

//============================================= destructor
dfem_diffusion::Solver::~Solver()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  MatDestroy(&A_);
}

//============================================= Initialize
void dfem_diffusion::Solver::Initialize()
{
  const std::string fname = "dfem_diffusion::Solver::Initialize";
  Chi::log.Log() << "\n"
                 << Chi::program_timer.GetTimeString() << " "
                 << TextName() << ": Initializing DFEM Diffusion solver ";

  //============================================= Get grid
  grid_ptr_ = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr_;
  if (grid_ptr_ == nullptr)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           " No grid defined.");

  Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= BIDs
  auto globl_unique_bndry_ids = grid.GetDomainUniqueBoundaryIDs();

  const auto& grid_boundary_id_map = grid_ptr_->GetBoundaryIDMap();
  for (uint64_t bndry_id : globl_unique_bndry_ids)
  {
    if (grid_boundary_id_map.count(bndry_id) == 0)
      throw std::logic_error(fname + ": Boundary id " +
                             std::to_string(bndry_id) + " does not have a name-assignment.");

    const auto& bndry_name = grid_boundary_id_map.at(bndry_id);
    if (boundary_preferences_.find(bndry_name) != boundary_preferences_.end())
    {
      BoundaryInfo bndry_info = boundary_preferences_.at(bndry_name);
      auto& bndry_vals = bndry_info.second;
      switch (bndry_info.first)
      {
        case BoundaryType::Reflecting:
        {
          boundaries_.insert(std::make_pair(
            bndry_id,Boundary{BoundaryType::Reflecting, {0., 0., 0.}}));
          Chi::log.Log() << "Boundary " << bndry_name << " set to reflecting.";
          break;
        }
        case BoundaryType::Dirichlet:
        {
          if (bndry_vals.empty()) bndry_vals.resize(1,0.0);
          boundaries_.insert(std::make_pair(
            bndry_id,Boundary{BoundaryType::Dirichlet, {bndry_vals[0], 0., 0.}}));
          Chi::log.Log() << "Boundary " << bndry_name << " set to dirichlet.";
          break;
        }
        case BoundaryType::Robin:
        {
          if (bndry_vals.size()!=3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Robin needs 3 values in bndry vals.");
          boundaries_.insert(std::make_pair(
            bndry_id,Boundary{BoundaryType::Robin, {bndry_vals[0],
                                                    bndry_vals[1],
                                                    bndry_vals[2]}}));
          Chi::log.Log() << "Boundary " << bndry_name
                         << " set to robin." << bndry_vals[0]<<","
                         << bndry_vals[1]<<","<<bndry_vals[2];
          break;
        }
        case BoundaryType::Vacuum:
        {
          boundaries_.insert(std::make_pair(
            bndry_id,Boundary{BoundaryType::Robin, {0.25, 0.5, 0.}}));
          Chi::log.Log() << "Boundary " << bndry_name << " set to vacuum.";
          break;
        }
        case BoundaryType::Neumann:
        {
          if (bndry_vals.size()!=3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Neumann needs 3 values in bndry vals.");
          boundaries_.insert(std::make_pair(
            bndry_id,Boundary{BoundaryType::Robin, {0., bndry_vals[0],
                                                    bndry_vals[1]}}));
          Chi::log.Log() << "Boundary " << bndry_name
                         << " set to neumann." << bndry_vals[0];
          break;
        }
      }//switch boundary type
    }
    else
    {
      boundaries_.insert(std::make_pair(
        bndry_id,Boundary{BoundaryType::Dirichlet, {0., 0., 0.}}));
      Chi::log.Log0Verbose1()
        << "No boundary preference found for boundary index " << bndry_name
        << "Dirichlet boundary added with zero boundary value.";
    }
  }//for bndry
  
  //============================================= Make SDM
  sdm_ptr_ = chi_math::spatial_discretization::PieceWiseLinearDiscontinuous::New(*grid_ptr_);
  const auto& sdm = *sdm_ptr_;
 
  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  num_local_dofs_ = sdm.GetNumLocalDOFs(OneDofPerNode);
  num_globl_dofs_ = sdm.GetNumGlobalDOFs(OneDofPerNode);

  Chi::log.Log() << "Num local DOFs: " << num_local_dofs_;
  Chi::log.Log() << "Num globl DOFs: " << num_globl_dofs_;

  //============================================= Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs_);
  const auto N = static_cast<int64_t>(num_globl_dofs_);

  A_ = chi_math::PETScUtils::CreateSquareMatrix(n, N);
  x_ = chi_math::PETScUtils::CreateVector(n, N);
  b_ = chi_math::PETScUtils::CreateVector(n, N);
 
  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag,nodal_nnz_off_diag, OneDofPerNode);
 
  chi_math::PETScUtils::InitMatrixSparsity(A_,
                                           nodal_nnz_in_diag,
                                           nodal_nnz_off_diag);

  if (field_functions_.empty())
  {
    std::string solver_name;
    if (not TextName().empty()) solver_name = TextName() + "-";

    std::string text_name = solver_name + "phi";

    using namespace chi_math;
    auto initial_field_function =
      std::make_shared<chi_physics::FieldFunctionGridBased>(
          text_name,                     //Text name
          sdm_ptr_,                       //Spatial Discretization
          Unknown(UnknownType::SCALAR)); //Unknown/Variable

    field_functions_.push_back(initial_field_function);
    Chi::field_function_stack.push_back(initial_field_function);
  }//if not ff set

}//end initialize

//========================================================== Execute
void dfem_diffusion::Solver::Execute()
{
  Chi::log.Log() << "\nExecuting DFEM IP Diffusion solver";

  const auto& grid = *grid_ptr_;
  const auto& sdm  = *sdm_ptr_;

  lua_State* L = Chi::console.GetConsoleState();

  //============================================= Assemble the system
  // is this needed?
  VecSet(b_, 0.0);

  Chi::log.Log() << "Assembling system: ";

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes   = cell_mapping.NumNodes();
    const auto   cc_nodes    = cell_mapping.GetNodeLocations();
    const auto  qp_data      = cell_mapping.MakeVolumetricQuadraturePointData();

    const auto imat  = cell.material_id_;
    MatDbl Acell(num_nodes, VecDbl(num_nodes, 0.0));
    VecDbl cell_rhs(num_nodes, 0.0);

    //==================================== Assemble volumetric terms
    for (size_t i=0; i<num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOF(cell, i);

      for (size_t j=0; j<num_nodes; ++j)
      {
        const int64_t jmap = sdm.MapDOF(cell, j);
        double entry_aij = 0.0;
        for (size_t qp : qp_data.QuadraturePointIndices())
        {
          entry_aij +=
            (
              CallLua_iXYZFunction(L, "D_coef",imat,qp_data.QPointXYZ(qp))  *
              qp_data.ShapeGrad(i, qp).Dot(qp_data.ShapeGrad(j, qp))
              +
              CallLua_iXYZFunction(L, "Sigma_a",imat,qp_data.QPointXYZ(qp))  *
              qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp)
            )
            *
            qp_data.JxW(qp);
        }//for qp
        MatSetValue(A_, imap, jmap, entry_aij, ADD_VALUES);
      }//for j
      double entry_rhs_i = 0.0;
      for (size_t qp : qp_data.QuadraturePointIndices())
        entry_rhs_i += CallLua_iXYZFunction(L,"Q_ext",imat,qp_data.QPointXYZ(qp))
          * qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      VecSetValue(b_, imap, entry_rhs_i, ADD_VALUES);
    }//for i

    //==================================== Assemble face terms
    const size_t num_faces = cell.faces_.size();
    for (size_t f = 0; f < num_faces; ++f) {
      const auto &face = cell.faces_[f];
      const auto &n_f = face.normal_;
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
      const auto fqp_data = cell_mapping.MakeSurfaceQuadraturePointData(f);

      const double hm = HPerpendicular(cell, f);

      typedef chi_mesh::MeshContinuum Grid;

      // interior face
      if (face.has_neighbor_)
      {
        const auto &adj_cell = grid.cells[face.neighbor_id_];
        const auto &adj_cell_mapping = sdm.GetCellMapping(adj_cell);
        const auto ac_nodes = adj_cell_mapping.GetNodeLocations();
        const size_t acf = Grid::MapCellFace(cell, adj_cell, f);
        const double hp_neigh = HPerpendicular(adj_cell, acf);

        const auto imat_neigh = adj_cell.material_id_;

        //========================= Compute Ckappa IP
        double Ckappa = 1.0;
        if (cell.Type() == chi_mesh::CellType::SLAB)
          Ckappa = 2.0;
        if (cell.Type() == chi_mesh::CellType::POLYGON)
          Ckappa = 2.0;
        if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
          Ckappa = 4.0;

        //========================= Assembly penalty terms
        for (size_t fi = 0; fi < num_face_nodes; ++fi) {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const int64_t imap = sdm.MapDOF(cell, i);

          for (size_t fj = 0; fj < num_face_nodes; ++fj) {
            const int jm = cell_mapping.MapFaceNode(f, fj);      //j-minus
            const int64_t jmmap = sdm.MapDOF(cell, jm);

            const int jp = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes,
                                           f, acf, fj);         //j-plus
            const int64_t jpmap = sdm.MapDOF(adj_cell, jp);

            double aij = 0.0;
            for (size_t qp: fqp_data.QuadraturePointIndices())
              aij += Ckappa *
                     (CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp)) / hm +
                      CallLua_iXYZFunction(L, "D_coef", imat_neigh, fqp_data.QPointXYZ(qp)) / hp_neigh) / 2.
                     *
                     fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) *
                     fqp_data.JxW(qp);

            MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A_, imap, jpmap, -aij, ADD_VALUES);
          }//for fj
        }//for fi

        //========================= Assemble gradient terms
        // For the following comments we use the notation:
        // Dk = 0.5* n dot nabla bk

        // {{D d_n b_i}}[[Phi]]
        // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-

        // loop over node of current cell (gradient of b_i)
        for (int i = 0; i < num_nodes; ++i) {
          const int64_t imap = sdm.MapDOF(cell, i);

          // loop over faces
          for (int fj = 0; fj < num_face_nodes; ++fj) {
            const int jm = cell_mapping.MapFaceNode(f, fj);      //j-minus
            const int64_t jmmap = sdm.MapDOF(cell, jm);
            const int jp = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes,
                                           f, acf, fj);         //j-plus
            const int64_t jpmap = sdm.MapDOF(adj_cell, jp);

            chi_mesh::Vector3 vec_aij;
            for (size_t qp: fqp_data.QuadraturePointIndices())
              vec_aij +=
                CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp)) *
                fqp_data.ShapeValue(jm, qp) * fqp_data.ShapeGrad(i, qp) *
                fqp_data.JxW(qp);
            const double aij = -0.5 * n_f.Dot(vec_aij);

            MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A_, imap, jpmap, -aij, ADD_VALUES);
          }//for fj
        }//for i

        // {{D d_n Phi}}[[b_i]]
        // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
        for (int fi = 0; fi < num_face_nodes; ++fi) {
          const int im = cell_mapping.MapFaceNode(f, fi);       //i-minus
          const int64_t immap = sdm.MapDOF(cell, im);

          const int ip = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes,
                                         f, acf, fi);            //i-plus
          const int64_t ipmap = sdm.MapDOF(adj_cell, ip);

          for (int j = 0; j < num_nodes; ++j) {
            const int64_t jmap = sdm.MapDOF(cell, j);

            chi_mesh::Vector3 vec_aij;
            for (size_t qp: fqp_data.QuadraturePointIndices())
              vec_aij +=
                CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp)) *
                fqp_data.ShapeValue(im, qp) * fqp_data.ShapeGrad(j, qp) *
                fqp_data.JxW(qp);
            const double aij = -0.5 * n_f.Dot(vec_aij);

            MatSetValue(A_, immap, jmap, aij, ADD_VALUES);
            MatSetValue(A_, ipmap, jmap, -aij, ADD_VALUES);
          }//for j
        }//for fi

      }//internal face
      else
      { // boundary face
        const auto &bndry = boundaries_[face.neighbor_id_];
        // Robin boundary
        if (bndry.type_ == BoundaryType::Robin)
        {
          const auto qp_face_data =
            cell_mapping.MakeSurfaceQuadraturePointData(f);

          const auto &aval = bndry.values_[0];
          const auto &bval = bndry.values_[1];
          const auto &fval = bndry.values_[2];

          Chi::log.Log0Verbose1() << "Boundary  set as Robin with a,b,f = ("
                         << aval << ","
                         << bval << ","
                         << fval << ") ";
          // true Robin when a!=0, otherwise, it is a Neumann:
          // Assert if b=0
          if (std::fabs(bval) < 1e-8)
            throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");

          for (size_t fi = 0; fi < num_face_nodes; fi++) {
            const uint i = cell_mapping.MapFaceNode(f, fi);
            const int64_t ir = sdm.MapDOF(cell, i);

            if (std::fabs(aval) >= 1.0e-12) {
              for (size_t fj = 0; fj < num_face_nodes; fj++) {
                const uint j = cell_mapping.MapFaceNode(f, fj);
                const int64_t jr = sdm.MapDOF(cell, j);

                double aij = 0.0;
                for (size_t qp: fqp_data.QuadraturePointIndices())
                  aij += fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(j, qp) *
                         fqp_data.JxW(qp);
                aij *= (aval / bval);

                MatSetValue(A_, ir, jr, aij, ADD_VALUES);
              }//for fj
            }//if a nonzero

            if (std::fabs(fval) >= 1.0e-12) {
              double rhs_val = 0.0;
              for (size_t qp: fqp_data.QuadraturePointIndices())
                rhs_val += fqp_data.ShapeValue(i, qp) * fqp_data.JxW(qp);
              rhs_val *= (fval / bval);

              VecSetValue(b_, ir, rhs_val, ADD_VALUES);
            }//if f nonzero
          }//for fi
        }//Robin BC
        else if (bndry.type_ == BoundaryType::Dirichlet) {
          const double bc_value = bndry.values_[0];
          //========================= Compute kappa
          double Ckappa = 2.0;
          if (cell.Type() == chi_mesh::CellType::SLAB)
            Ckappa = 4.0; // fmax(4.0*Dg/hm,0.25);
          if (cell.Type() == chi_mesh::CellType::POLYGON)
            Ckappa = 4.0;
          if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
            Ckappa = 8.0;

          //========================= Assembly penalty terms
          for (size_t fi = 0; fi < num_face_nodes; ++fi) {
            const uint i = cell_mapping.MapFaceNode(f, fi);
            const int64_t imap = sdm.MapDOF(cell, i);

            for (size_t fj = 0; fj < num_face_nodes; ++fj) {
              const uint jm = cell_mapping.MapFaceNode(f, fj);
              const int64_t jmmap = sdm.MapDOF(cell, jm);

              double aij = 0.0;
              for (size_t qp: fqp_data.QuadraturePointIndices())
                aij += Ckappa *
                       CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp)) / hm *
                       fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) *
                       fqp_data.JxW(qp);
              double aij_bc_value = aij * bc_value;

              MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
              VecSetValue(b_, imap, aij_bc_value, ADD_VALUES);
            }//for fj
          }//for fi

          //========================= Assemble gradient terms
          // For the following comments we use the notation:
          // Dk = 0.5* n dot nabla bk

          // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
          for (size_t i = 0; i < num_nodes; i++) {
            const int64_t imap = sdm.MapDOF(cell, i);

            for (size_t j = 0; j < num_nodes; j++) {
              const int64_t jmap = sdm.MapDOF(cell, j);

              chi_mesh::Vector3 vec_aij;
              for (size_t qp: fqp_data.QuadraturePointIndices())
                vec_aij +=
                  (fqp_data.ShapeValue(j, qp) * fqp_data.ShapeGrad(i, qp) +
                   fqp_data.ShapeValue(i, qp) * fqp_data.ShapeGrad(j, qp)) *
                  fqp_data.JxW(qp) *
                  CallLua_iXYZFunction(L, "D_coef", imat, fqp_data.QPointXYZ(qp));

              const double aij = -n_f.Dot(vec_aij);
              double aij_bc_value = aij * bc_value;

              MatSetValue(A_, imap, jmap, aij, ADD_VALUES);
              VecSetValue(b_, imap, aij_bc_value, ADD_VALUES);
            }//for fj
          }//for i
        }//Dirichlet BC
        else {

        }//else BC
      }//boundary face
    }//for face f
  }//for cell

  Chi::log.Log() << "Global assembly";
 
  MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);

//  MatView(A, PETSC_VIEWER_STDERR_WORLD);
//
//  PetscViewer viewer;
//  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"A.m",&viewer);
//  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
//  MatView(A,viewer);
//  PetscViewerPopFormat(viewer);
//  PetscViewerDestroy(&viewer);

  Chi::log.Log() << "Done global assembly";

  //============================================= Create Krylov Solver
  Chi::log.Log() << "Solving: ";
  auto petsc_solver =
    chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
        A_,               //Matrix
        TextName(),      //Solver name
        KSPCG,           //Solver type
        PCGAMG,          //Preconditioner type
        basic_options_("residual_tolerance").FloatValue(),  //Relative residual tolerance
        basic_options_("max_iters").IntegerValue()          //Max iterations
    );

  //============================================= Solve
  KSPSolve(petsc_solver.ksp, b_, x_);

  Chi::log.Log() << "Done solving";

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  sdm.LocalizePETScVector(x_, field_, OneDofPerNode);

  field_functions_.front()->UpdateFieldVector(field_);
}