// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion_pwlc_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_structs.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"

namespace opensn
{

DiffusionPWLCSolver::DiffusionPWLCSolver(std::string text_name,
                                         const opensn::SpatialDiscretization& sdm,
                                         const UnknownManager& uk_man,
                                         std::map<uint64_t, BoundaryCondition> bcs,
                                         MatID2XSMap map_mat_id_2_xs,
                                         const std::vector<UnitCellMatrices>& unit_cell_matrices,
                                         const bool suppress_bcs,
                                         const bool verbose)
  : DiffusionSolver(std::move(text_name),
                    sdm,
                    uk_man,
                    std::move(bcs),
                    std::move(map_mat_id_2_xs),
                    unit_cell_matrices,
                    suppress_bcs,
                    true,
                    verbose)
{
  if (sdm_.Type() != SpatialDiscretizationType::PIECEWISE_LINEAR_CONTINUOUS)
    throw std::logic_error("acceleration::DiffusionPWLCSolver can only be used with PWLC.");
}

void
DiffusionPWLCSolver::AssembleAand_b(const std::vector<double>& q_vector)
{
  const size_t num_local_dofs = sdm_.GetNumLocalAndGhostDOFs(uk_man_);
  OpenSnInvalidArgumentIf(q_vector.size() != num_local_dofs,
                          std::string("q_vector size mismatch. ") +
                            std::to_string(q_vector.size()) + " vs " +
                            std::to_string(num_local_dofs));

  const std::string fname = "acceleration::DiffusionMIPSolver::"
                            "AssembleAand_b";
  if (A_ == nullptr or rhs_ == nullptr or ksp_ == nullptr)
    throw std::logic_error(fname + ": Some or all PETSc elements are null. "
                                   "Check that Initialize has been called.");
  if (options.verbose)
    log.Log() << program_timer.GetTimeString() << " Starting assembly";

  const size_t num_groups = uk_man_.unknowns.front().num_components;

  VecSet(rhs_, 0.0);
  for (const auto& cell : grid_.local_cells)
  {
    const size_t num_faces = cell.faces.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto cc_nodes = cell_mapping.GetNodeLocations();
    const auto& unit_cell_matrices = unit_cell_matrices_[cell.local_id];

    const auto& intV_gradshapeI_gradshapeJ = unit_cell_matrices.intV_gradshapeI_gradshapeJ;
    const auto& intV_shapeI_shapeJ = unit_cell_matrices.intV_shapeI_shapeJ;

    const auto& xs = mat_id_2_xs_map_.at(cell.material_id);

    // Mark dirichlet nodes
    std::vector<std::pair<bool, double>> node_is_dirichlet(num_nodes, {false, 0.0});
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      if (not face.has_neighbor and not suppress_bcs_)
      {
        BoundaryCondition bc;
        if (bcs_.count(face.neighbor_id) > 0)
          bc = bcs_.at(face.neighbor_id);

        if (bc.type != BCType::DIRICHLET)
          continue;

        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
          node_is_dirichlet[cell_mapping.MapFaceNode(f, fi)] = {true, bc.values[0]};
      }
    }

    for (size_t g = 0; g < num_groups; ++g)
    {
      // Get coefficient and nodal src
      const double Dg = xs.Dg[g];
      const double sigr_g = xs.sigR[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; j++)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      // Assemble continuous terms
      for (size_t i = 0; i < num_nodes; i++)
      {
        if (node_is_dirichlet[i].first)
          continue;
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0;
        for (size_t j = 0; j < num_nodes; j++)
        {
          const int64_t jmap = sdm_.MapDOF(cell, j, uk_man_, 0, g);

          const double entry_aij =
            Dg * intV_gradshapeI_gradshapeJ[i][j] + sigr_g * intV_shapeI_shapeJ[i][j];

          if (not node_is_dirichlet[j].first)
            MatSetValue(A_, imap, jmap, entry_aij, ADD_VALUES);
          else
          {
            const double bcvalue = node_is_dirichlet[j].second;
            VecSetValue(rhs_, imap, -entry_aij * bcvalue, ADD_VALUES);
          }

          entry_rhs_i += intV_shapeI_shapeJ[i][j] * qg[j];
        } // for j

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      // Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& intS_shapeI_shapeJ = unit_cell_matrices.intS_shapeI_shapeJ[f];
        const auto& intS_shapeI = unit_cell_matrices.intS_shapeI[f];

        if (not face.has_neighbor and not suppress_bcs_)
        {
          BoundaryCondition bc;
          if (bcs_.count(face.neighbor_id) > 0)
            bc = bcs_.at(face.neighbor_id);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              // MatSetValue(A_, imap, imap, intV_shapeI[i], ADD_VALUES);
              // VecSetValue(rhs_, imap, bc_value * intV_shapeI[i], ADD_VALUES);
              MatSetValue(A_, imap, imap, 1.0, ADD_VALUES);
              VecSetValue(rhs_, imap, bc_value, ADD_VALUES);
            } // for fi

          } // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double aval = bc.values[0];
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12)
              continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; fi++)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(aval) >= 1.0e-12)
              {
                for (size_t fj = 0; fj < num_face_nodes; fj++)
                {
                  const int j = cell_mapping.MapFaceNode(f, fj);
                  const int64_t jr = sdm_.MapDOF(cell, j, uk_man_, 0, g);

                  const double aij = (aval / bval) * intS_shapeI_shapeJ[i][j];

                  MatSetValue(A_, ir, jr, aij, ADD_VALUES);
                } // for fj
              }   // if a nonzero

              if (std::fabs(fval) >= 1.0e-12)
              {
                const double rhs_val = (fval / bval) * intS_shapeI[i];

                VecSetValue(rhs_, ir, rhs_val, ADD_VALUES);
              } // if f nonzero
            }   // for fi
          }     // Robin BC
        }       // boundary face
      }         // for face
    }           // for g
  }             // for cell

  MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(rhs_);
  VecAssemblyEnd(rhs_);

  if (options.verbose)
  {
    MatInfo info;
    MatGetInfo(A_, MAT_GLOBAL_SUM, &info);

    log.Log() << "Number of mallocs used = " << info.mallocs
              << "\nNumber of non-zeros allocated = " << info.nz_allocated
              << "\nNumber of non-zeros used = " << info.nz_used
              << "\nNumber of unneeded non-zeros = " << info.nz_unneeded;
  }

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(A_, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  KSPSetOperators(ksp_, A_, A_);

  if (options.verbose)
    log.Log() << program_timer.GetTimeString() << " Assembly completed";

  PC pc;
  KSPGetPC(ksp_, &pc);
  PCSetUp(pc);

  KSPSetUp(ksp_);
}

void
DiffusionPWLCSolver::Assemble_b(const std::vector<double>& q_vector)
{
  const size_t num_local_dofs = sdm_.GetNumLocalAndGhostDOFs(uk_man_);
  OpenSnInvalidArgumentIf(q_vector.size() != num_local_dofs,
                          std::string("q_vector size mismatch. ") +
                            std::to_string(q_vector.size()) + " vs " +
                            std::to_string(num_local_dofs));
  const std::string fname = "acceleration::DiffusionMIPSolver::"
                            "Assemble_b";
  if (A_ == nullptr or rhs_ == nullptr or ksp_ == nullptr)
    throw std::logic_error(fname + ": Some or all PETSc elements are null. "
                                   "Check that Initialize has been called.");
  if (options.verbose)
    log.Log() << program_timer.GetTimeString() << " Starting assembly";

  const size_t num_groups = uk_man_.unknowns.front().num_components;

  VecSet(rhs_, 0.0);
  for (const auto& cell : grid_.local_cells)
  {
    const size_t num_faces = cell.faces.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto cc_nodes = cell_mapping.GetNodeLocations();
    const auto& unit_cell_matrices = unit_cell_matrices_[cell.local_id];

    const auto& intV_gradshapeI_gradshapeJ = unit_cell_matrices.intV_gradshapeI_gradshapeJ;
    const auto& intV_shapeI_shapeJ = unit_cell_matrices.intV_shapeI_shapeJ;
    const auto& intV_shapeI = unit_cell_matrices.intV_shapeI;

    const auto& xs = mat_id_2_xs_map_.at(cell.material_id);

    // Mark dirichlet nodes
    std::vector<std::pair<bool, double>> node_is_dirichlet(num_nodes, {false, 0.0});
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      if (not face.has_neighbor and suppress_bcs_)
      {
        BoundaryCondition bc;
        if (bcs_.count(face.neighbor_id) > 0)
          bc = bcs_.at(face.neighbor_id);

        if (bc.type != BCType::DIRICHLET)
          continue;

        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
          node_is_dirichlet[cell_mapping.MapFaceNode(f, fi)] = {true, bc.values[0]};
      }
    }

    for (size_t g = 0; g < num_groups; ++g)
    {
      // Get coefficient and nodal src
      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; j++)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      // Assemble continuous terms
      const double Dg = xs.Dg[g];
      const double sigr_g = xs.sigR[g];

      for (size_t i = 0; i < num_nodes; i++)
      {
        if (node_is_dirichlet[i].first)
          continue;
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0;
        for (size_t j = 0; j < num_nodes; j++)
        {
          if (node_is_dirichlet[j].first)
          {
            const double entry_aij =
              Dg * intV_gradshapeI_gradshapeJ[i][j] + sigr_g * intV_shapeI_shapeJ[i][j];

            const double bcvalue = node_is_dirichlet[j].second;
            VecSetValue(rhs_, imap, -entry_aij * bcvalue, ADD_VALUES);
          }

          entry_rhs_i += intV_shapeI_shapeJ[i][j] * qg[j];
        } // for j

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      // Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& intS_shapeI = unit_cell_matrices.intS_shapeI[f];

        if (not face.has_neighbor and suppress_bcs_)
        {
          BoundaryCondition bc;
          if (bcs_.count(face.neighbor_id) > 0)
            bc = bcs_.at(face.neighbor_id);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            // Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              // VecSetValue(rhs_, imap, bc_value * intV_shapeI[i], ADD_VALUES);
              VecSetValue(rhs_, imap, bc_value, ADD_VALUES);
            } // for fi

          } // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12)
              continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; fi++)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(fval) >= 1.0e-12)
              {
                const double rhs_val = (fval / bval) * intS_shapeI[i];

                VecSetValue(rhs_, ir, rhs_val, ADD_VALUES);
              } // if f nonzero
            }   // for fi
          }     // Robin BC
        }       // boundary face
      }         // for face
    }           // for g
  }             // for cell

  VecAssemblyBegin(rhs_);
  VecAssemblyEnd(rhs_);

  if (options.verbose)
    log.Log() << program_timer.GetTimeString() << " Assembly completed";
}

void
DiffusionPWLCSolver::Assemble_b(Vec petsc_q_vector)
{
  const std::string fname = "acceleration::DiffusionMIPSolver::"
                            "Assemble_b";
  if (A_ == nullptr or rhs_ == nullptr or ksp_ == nullptr)
    throw std::logic_error(fname + ": Some or all PETSc elements are null. "
                                   "Check that Initialize has been called.");
  if (options.verbose)
    log.Log() << program_timer.GetTimeString() << " Starting assembly";

  const size_t num_groups = uk_man_.unknowns.front().num_components;

  const double* q_vector;
  VecGetArrayRead(petsc_q_vector, &q_vector);

  VecSet(rhs_, 0.0);
  for (const auto& cell : grid_.local_cells)
  {
    const size_t num_faces = cell.faces.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto cc_nodes = cell_mapping.GetNodeLocations();
    const auto& unit_cell_matrices = unit_cell_matrices_[cell.local_id];

    const auto& intV_shapeI_shapeJ = unit_cell_matrices.intV_shapeI_shapeJ;
    const auto& intV_shapeI = unit_cell_matrices.intV_shapeI;

    // Mark dirichlet nodes
    std::vector<bool> node_is_dirichlet(num_nodes, false);
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      if (not face.has_neighbor and not suppress_bcs_)
      {
        BoundaryCondition bc;
        if (bcs_.count(face.neighbor_id) > 0)
          bc = bcs_.at(face.neighbor_id);

        if (bc.type != BCType::DIRICHLET)
          continue;

        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
          node_is_dirichlet[cell_mapping.MapFaceNode(f, fi)] = true;
      }
    }

    for (size_t g = 0; g < num_groups; ++g)
    {
      // Get coefficient and nodal src
      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; j++)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      // Assemble continuous terms
      for (size_t i = 0; i < num_nodes; i++)
      {
        if (node_is_dirichlet[i])
          continue;
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0;
        for (size_t j = 0; j < num_nodes; j++)
          entry_rhs_i += intV_shapeI_shapeJ[i][j] * qg[j];

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      // Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& intS_shapeI = unit_cell_matrices.intS_shapeI[f];

        if (not face.has_neighbor and not suppress_bcs_)
        {
          BoundaryCondition bc;
          if (bcs_.count(face.neighbor_id) > 0)
            bc = bcs_.at(face.neighbor_id);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            // Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              VecSetValue(rhs_, imap, bc_value * intV_shapeI[i], ADD_VALUES);
            } // for fi

          } // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12)
              continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; fi++)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(fval) >= 1.0e-12)
              {
                const double rhs_val = (fval / bval) * intS_shapeI[i];

                VecSetValue(rhs_, ir, rhs_val, ADD_VALUES);
              } // if f nonzero
            }   // for fi
          }     // Robin BC
        }       // boundary face
      }         // for face
    }           // for g
  }             // for cell

  VecRestoreArrayRead(petsc_q_vector, &q_vector);

  VecAssemblyBegin(rhs_);
  VecAssemblyEnd(rhs_);

  if (options.verbose)
    log.Log() << program_timer.GetTimeString() << " Assembly completed";
}

} // namespace opensn
