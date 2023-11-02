#include "diffusion_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

// ###################################################################
/**Assembles PWLC matrix for general cells.*/
void chi_diffusion::Solver::CFEM_Assemble_A_and_b(chi_mesh::Cell& cell,
                                                  int group)
{
  auto pwl_sdm = std::static_pointer_cast<
    chi_math::spatial_discretization::PieceWiseLinearContinuous>(this->discretization_);
  const auto& fe_intgrl_values = unit_integrals_.at(cell.global_id_);

  size_t num_nodes = fe_intgrl_values.NumNodes();

  //======================================== Process material id
  std::vector<double> D(num_nodes, 1.0);
  std::vector<double> q(num_nodes, 1.0);
  std::vector<double> siga(num_nodes, 0.0);

  GetMaterialProperties(cell, num_nodes, D, q, siga, group);

  //======================================== Init cell matrix info
  typedef std::vector<double> Row;
  typedef std::vector<Row> Matrix;

  Matrix              cell_matrix;
  std::vector<double> cell_rhs;

  cell_matrix.resize(num_nodes, Row(num_nodes, 0.0));
  cell_rhs.resize(num_nodes, 0.0);

  std::vector<int64_t> dof_global_row_ind(num_nodes, -1);
  std::vector<int64_t> dof_global_col_ind(num_nodes, -1);

  //========================================= Loop over DOFs
  for (int i=0; i<num_nodes; i++)
  {
    dof_global_row_ind[i] = pwl_sdm->MapDOF(cell,i);

    for (int j=0; j<num_nodes; j++)
    {
      double mat_entry =
        D[j]* fe_intgrl_values.IntV_gradShapeI_gradShapeJ(i, j) +
        siga[j]* fe_intgrl_values.IntV_shapeI_shapeJ(i, j);

      cell_matrix[i][j] = mat_entry;
    }//for j

    //====================== Develop RHS entry
    cell_rhs[i] = q[i]* fe_intgrl_values.IntV_shapeI(i);
  }//for i
  dof_global_col_ind = dof_global_row_ind;

  //======================================== Apply Dirichlet,Vacuum, Neumann and
  //                                         Robin BCs
  // Dirichlets are just collected
  std::vector<int>    dirichlet_count(num_nodes, 0);
  std::vector<double> dirichlet_value(num_nodes, 0.0);
  for (int f=0; f<cell.faces_.size(); f++)
  {
    if (not cell.faces_[f].has_neighbor_)
    {
      uint64_t ir_boundary_index = cell.faces_[f].neighbor_id_;
      auto ir_boundary_type  = boundaries_.at(ir_boundary_index)->type_;

      if (ir_boundary_type == BoundaryType::Dirichlet)
      {
        auto& dirichlet_bndry =
          (chi_diffusion::BoundaryDirichlet&)*boundaries_.at(ir_boundary_index);

        int num_face_dofs = cell.faces_[f].vertex_ids_.size();
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i  = fe_intgrl_values.FaceDofMapping(f,fi);
          dirichlet_count[i] += 1;
          dirichlet_value[i] += dirichlet_bndry.boundary_value;
        }
      }

      if (ir_boundary_type == BoundaryType::Robin)
      {
        auto& robin_bndry =
          (chi_diffusion::BoundaryRobin&)*boundaries_.at(ir_boundary_index);

        std::cout << robin_bndry.a << " " << robin_bndry.b << " " << robin_bndry.f << std::endl;

        int num_face_dofs = cell.faces_[f].vertex_ids_.size();
        for (int fi=0; fi<num_face_dofs; fi++)
        {
          int i  = fe_intgrl_values.FaceDofMapping(f,fi);

          for (int fj=0; fj<num_face_dofs; fj++)
          {
            int j  = fe_intgrl_values.FaceDofMapping(f,fj);

            double aij = robin_bndry.a* fe_intgrl_values.IntS_shapeI_shapeJ(f, i, j);
            aij /= robin_bndry.b;

            cell_matrix[i][j] += aij;
          }//for fj

          double aii = robin_bndry.f* fe_intgrl_values.IntS_shapeI(f, i);
          aii /= robin_bndry.b;

          cell_matrix[i][i] += aii;
        }//for fi
      }//if robin

    }//if boundary

  }//for face

  //======================================== Apply dirichlet BCs
  //Compute average dirichlet value
  for (int i=0; i<num_nodes; ++i)
    dirichlet_value[i] /= (dirichlet_count[i] > 0)? dirichlet_count[i] : 1;

  for (int i=0; i<num_nodes; ++i)
  {
    if (dirichlet_count[i] > 0)
    {
      cell_matrix[i].clear();
      cell_matrix[i] = std::vector<double>(num_nodes, 0.0);
      cell_matrix[i][i] = 1.0;
      int ir = dof_global_col_ind[i];
      MatSetValue(A_, ir, ir, 1.0, ADD_VALUES);
      dof_global_col_ind[i] = -1;
      cell_rhs[i] = dirichlet_value[i];
    }
    else
    {
      for (int j=0; j<num_nodes; ++j)
      {
        if (dirichlet_count[j] > 0)
        {
          cell_rhs[i] -= cell_matrix[i][j]*dirichlet_value[j];
          cell_matrix[i][j] = 0.0;
        }
      }
    }
  }

  //======================================== Make contiguous copy of matrix
  std::vector<double> cell_matrix_cont(num_nodes * num_nodes, 0.0);
  int n = 0;
  for (int i=0; i<num_nodes; ++i)
    for (int j=0; j<num_nodes; ++j)
      cell_matrix_cont[n++] = cell_matrix[i][j];

  //======================================== Add to global
  MatSetValues(A_,
               num_nodes, dof_global_row_ind.data(),
               num_nodes, dof_global_col_ind.data(),
               cell_matrix_cont.data(), ADD_VALUES);

  VecSetValues(b_,
               num_nodes, dof_global_row_ind.data(),
               cell_rhs.data(), ADD_VALUES);

  VecSetValues(x_,
               num_nodes,
               dof_global_row_ind.data(),
               dirichlet_value.data(), INSERT_VALUES);

}
