// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_structs.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/functions/scalar_spatial_function.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include <utility>

namespace opensn
{

void
DiffusionMIPSolver::SetSourceFunction(std::shared_ptr<ScalarSpatialFunction> function)
{
  source_function_ = function;
}

void
DiffusionMIPSolver::SetReferenceSolutionFunction(std::shared_ptr<ScalarSpatialFunction> function)
{
  ref_solution_function_ = function;
}

DiffusionMIPSolver::DiffusionMIPSolver(std::string name,
                                       const opensn::SpatialDiscretization& sdm,
                                       const UnknownManager& uk_man,
                                       std::map<uint64_t, BoundaryCondition> bcs,
                                       MatID2XSMap map_mat_id_2_xs,
                                       const std::vector<UnitCellMatrices>& unit_cell_matrices,
                                       const bool suppress_bcs,
                                       const bool verbose)
  : DiffusionSolver(std::move(name),
                    sdm,
                    uk_man,
                    std::move(bcs),
                    std::move(map_mat_id_2_xs),
                    unit_cell_matrices,
                    false,
                    suppress_bcs,
                    verbose)
{
  if (sdm_.Type() != SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::logic_error("acceleration::DiffusionMIPSolver can only be used with PWLD.");
}

void
DiffusionMIPSolver::AssembleAand_b_wQpoints(const std::vector<double>& q_vector)
{
  const std::string fname = "acceleration::DiffusionMIPSolver::"
                            "AssembleAand_b_wQpoints";
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
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    const auto& xs = mat_id_2_xs_map_.at(cell.material_id);

    // For component/group
    for (size_t g = 0; g < num_groups; ++g)
    {
      // Get coefficient and nodal src
      const double Dg = xs.Dg[g];
      const double sigr_g = xs.sigR[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; ++j)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      // Assemble continuous terms
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0; // entry may accumulate over j
        for (size_t j = 0; j < num_nodes; ++j)
        {
          const int64_t jmap = sdm_.MapDOF(cell, j, uk_man_, 0, g);
          double entry_aij = 0.0;
          for (size_t qp : fe_vol_data.QuadraturePointIndices())
          {
            entry_aij += Dg * fe_vol_data.ShapeGrad(i, qp).Dot(fe_vol_data.ShapeGrad(j, qp)) *
                         fe_vol_data.JxW(qp);

            entry_aij += sigr_g * fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeValue(j, qp) *
                         fe_vol_data.JxW(qp);

            if (not source_function_)
              entry_rhs_i += fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeValue(j, qp) *
                             fe_vol_data.JxW(qp) * qg[j];
          } // for qp
          MatSetValue(A_, imap, jmap, entry_aij, ADD_VALUES);
        } // for j

        if (source_function_)
        {
          for (size_t qp : fe_vol_data.QuadraturePointIndices())
            entry_rhs_i += source_function_->Evaluate(fe_vol_data.QPointXYZ(qp)) *
                           fe_vol_data.ShapeValue(i, qp) * fe_vol_data.JxW(qp);
        }

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      // Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const auto& n_f = face.normal;
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);

        const double hm = HPerpendicular(cell, f);

        if (face.has_neighbor)
        {
          const auto& adj_cell = grid_.cells[face.neighbor_id];
          const auto& adj_cell_mapping = sdm_.GetCellMapping(adj_cell);
          const auto ac_nodes = adj_cell_mapping.GetNodeLocations();
          const size_t acf = MeshContinuum::MapCellFace(cell, adj_cell, f);
          const double hp = HPerpendicular(adj_cell, acf);

          const auto& adj_xs = mat_id_2_xs_map_.at(adj_cell.material_id);
          const double adj_Dg = adj_xs.Dg[g];

          // Compute kappa
          double kappa = 1.0;
          if (cell.Type() == CellType::SLAB)
            kappa = fmax(options.penalty_factor * (adj_Dg / hp + Dg / hm) * 0.5, 0.25);
          if (cell.Type() == CellType::POLYGON)
            kappa = fmax(options.penalty_factor * (adj_Dg / hp + Dg / hm) * 0.5, 0.25);
          if (cell.Type() == CellType::POLYHEDRON)
            kappa = fmax(options.penalty_factor * 2.0 * (adj_Dg / hp + Dg / hm) * 0.5, 0.25);

          // Assembly penalty terms
          for (size_t fi = 0; fi < num_face_nodes; ++fi)
          {
            const int i = cell_mapping.MapFaceNode(f, fi);
            const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

            for (size_t fj = 0; fj < num_face_nodes; ++fj)
            {
              const int jm = cell_mapping.MapFaceNode(f, fj); // j-minus
              const int jp =
                MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fj); // j-plus
              const int64_t jmmap = sdm_.MapDOF(cell, jm, uk_man_, 0, g);
              const int64_t jpmap = sdm_.MapDOF(adj_cell, jp, uk_man_, 0, g);

              double aij = 0.0;
              for (size_t qp : fe_srf_data.QuadraturePointIndices())
                aij += kappa * fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(jm, qp) *
                       fe_srf_data.JxW(qp);

              MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
              MatSetValue(A_, imap, jpmap, -aij, ADD_VALUES);
            } // for fj
          }   // for fi

          // Assemble gradient terms
          // For the following comments we use the notation:
          // Dk = 0.5* n dot nabla bk

          // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
          for (int i = 0; i < num_nodes; ++i)
          {
            const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

            for (int fj = 0; fj < num_face_nodes; ++fj)
            {
              const int jm = cell_mapping.MapFaceNode(f, fj); // j-minus
              const int jp =
                MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fj); // j-plus
              const int64_t jmmap = sdm_.MapDOF(cell, jm, uk_man_, 0, g);
              const int64_t jpmap = sdm_.MapDOF(adj_cell, jp, uk_man_, 0, g);

              Vector3 vec_aij;
              for (size_t qp : fe_srf_data.QuadraturePointIndices())
                vec_aij += fe_srf_data.ShapeValue(jm, qp) * fe_srf_data.ShapeGrad(i, qp) *
                           fe_srf_data.JxW(qp);
              const double aij = -0.5 * Dg * n_f.Dot(vec_aij);

              MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
              MatSetValue(A_, imap, jpmap, -aij, ADD_VALUES);
            } // for fj
          }   // for i

          // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
          for (int fi = 0; fi < num_face_nodes; ++fi)
          {
            const int im = cell_mapping.MapFaceNode(f, fi); // i-minus
            const int ip =
              MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fi); // i-plus
            const int64_t immap = sdm_.MapDOF(cell, im, uk_man_, 0, g);
            const int64_t ipmap = sdm_.MapDOF(adj_cell, ip, uk_man_, 0, g);

            for (int j = 0; j < num_nodes; ++j)
            {
              const int64_t jmap = sdm_.MapDOF(cell, j, uk_man_, 0, g);

              Vector3 vec_aij;
              for (size_t qp : fe_srf_data.QuadraturePointIndices())
                vec_aij += fe_srf_data.ShapeValue(im, qp) * fe_srf_data.ShapeGrad(j, qp) *
                           fe_srf_data.JxW(qp);
              const double aij = -0.5 * Dg * n_f.Dot(vec_aij);

              MatSetValue(A_, immap, jmap, aij, ADD_VALUES);
              MatSetValue(A_, ipmap, jmap, -aij, ADD_VALUES);
            } // for j
          }   // for fi

        } // internal face
        else
        {
          BoundaryCondition bc;
          if (bcs_.count(face.neighbor_id) > 0)
            bc = bcs_.at(face.neighbor_id);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            // Compute kappa
            double kappa = 1.0;
            if (cell.Type() == CellType::SLAB)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYGON)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYHEDRON)
              kappa = fmax(options.penalty_factor * 2.0 * Dg / hm, 0.25);

            // Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const int jm = cell_mapping.MapFaceNode(f, fj);
                const int64_t jmmap = sdm_.MapDOF(cell, jm, uk_man_, 0, g);

                double aij = 0.0;
                for (size_t qp : fe_srf_data.QuadraturePointIndices())
                  aij += kappa * fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(jm, qp) *
                         fe_srf_data.JxW(qp);
                double aij_bc_value = aij * bc_value;

                if (ref_solution_function_)
                {
                  aij_bc_value = 0.0;
                  for (size_t qp : fe_srf_data.QuadraturePointIndices())
                    aij_bc_value += kappa *
                                    ref_solution_function_->Evaluate(fe_srf_data.QPointXYZ(qp)) *
                                    fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(jm, qp) *
                                    fe_srf_data.JxW(qp);
                }

                MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for fi

            // Assemble gradient terms
            // For the following comments we use the notation:
            // Dk = n dot nabla bk

            // D* n dot (b_j^+ - b_j^-)*nabla b_i^-
            for (size_t i = 0; i < num_nodes; ++i)
            {
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t j = 0; j < num_nodes; ++j)
              {
                const int64_t jmap = sdm_.MapDOF(cell, j, uk_man_, 0, g);

                Vector3 vec_aij;
                for (size_t qp : fe_srf_data.QuadraturePointIndices())
                  vec_aij += fe_srf_data.ShapeValue(j, qp) * fe_srf_data.ShapeGrad(i, qp) *
                               fe_srf_data.JxW(qp) +
                             fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeGrad(j, qp) *
                               fe_srf_data.JxW(qp);
                const double aij = -Dg * n_f.Dot(vec_aij);

                double aij_bc_value = aij * bc_value;

                if (ref_solution_function_)
                {
                  Vector3 vec_aij_mms;
                  for (size_t qp : fe_srf_data.QuadraturePointIndices())
                    vec_aij_mms += ref_solution_function_->Evaluate(fe_srf_data.QPointXYZ(qp)) *
                                   (fe_srf_data.ShapeValue(j, qp) * fe_srf_data.ShapeGrad(i, qp) *
                                      fe_srf_data.JxW(qp) +
                                    fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeGrad(j, qp) *
                                      fe_srf_data.JxW(qp));
                  aij_bc_value = -Dg * n_f.Dot(vec_aij_mms);
                }

                MatSetValue(A_, imap, jmap, aij, ADD_VALUES);
                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for i
          }     // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double aval = bc.values[0];
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12)
              continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(aval) >= 1.0e-12)
              {
                for (size_t fj = 0; fj < num_face_nodes; ++fj)
                {
                  const int j = cell_mapping.MapFaceNode(f, fj);
                  const int64_t jr = sdm_.MapDOF(cell, j, uk_man_, 0, g);

                  double aij = 0.0;
                  for (size_t qp : fe_srf_data.QuadraturePointIndices())
                    aij += fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(j, qp) *
                           fe_srf_data.JxW(qp);
                  aij *= (aval / bval);

                  MatSetValue(A_, ir, jr, aij, ADD_VALUES);
                } // for fj
              }   // if a nonzero

              if (std::fabs(fval) >= 1.0e-12)
              {
                double rhs_val = 0.0;
                for (size_t qp : fe_srf_data.QuadraturePointIndices())
                  rhs_val += fe_srf_data.ShapeValue(i, qp) * fe_srf_data.JxW(qp);
                rhs_val *= (fval / bval);

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
DiffusionMIPSolver::Assemble_b_wQpoints(const std::vector<double>& q_vector)
{
  const std::string fname = "acceleration::DiffusionMIPSolver::"
                            "AssembleAand_b_wQpoints";
  if (A_ == nullptr or rhs_ == nullptr or ksp_ == nullptr)
    throw std::logic_error(fname + ": Some or all PETSc elements are null. "
                                   "Check that Initialize has been called.");
  if (options.verbose)
    log.Log() << program_timer.GetTimeString() << " Starting assembly";

  VecSet(rhs_, 0.0);

  for (const auto& cell : grid_.local_cells)
  {
    const size_t num_faces = cell.faces.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();
    const size_t num_groups = uk_man_.unknowns.front().num_components;

    const auto& xs = mat_id_2_xs_map_.at(cell.material_id);

    // For component/group
    for (size_t g = 0; g < num_groups; ++g)
    {
      // Get coefficient and nodal src
      const double Dg = xs.Dg[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; ++j)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      // Assemble continuous terms
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0; // entry may accumulate over j
        if (not source_function_)
          for (size_t j = 0; j < num_nodes; ++j)
          {
            for (size_t qp : fe_vol_data.QuadraturePointIndices())
            {
              entry_rhs_i += fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeValue(j, qp) *
                             fe_vol_data.JxW(qp) * qg[j];
            } // for qp
          }   // for j
        else
        {
          for (size_t qp : fe_vol_data.QuadraturePointIndices())
            entry_rhs_i += source_function_->Evaluate(fe_vol_data.QPointXYZ(qp)) *
                           fe_vol_data.ShapeValue(i, qp) * fe_vol_data.JxW(qp);
        }

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      // Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const auto& n_f = face.normal;
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);

        const double hm = HPerpendicular(cell, f);

        if (not face.has_neighbor and not suppress_bcs_)
        {
          BoundaryCondition bc;
          if (bcs_.count(face.neighbor_id) > 0)
            bc = bcs_.at(face.neighbor_id);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            // Compute kappa
            double kappa = 1.0;
            if (cell.Type() == CellType::SLAB)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYGON)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYHEDRON)
              kappa = fmax(options.penalty_factor * 2.0 * Dg / hm, 0.25);

            // Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const int jm = cell_mapping.MapFaceNode(f, fj);

                double aij = 0.0;
                for (size_t qp : fe_srf_data.QuadraturePointIndices())
                  aij += kappa * fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(jm, qp) *
                         fe_srf_data.JxW(qp);
                double aij_bc_value = aij * bc_value;

                if (ref_solution_function_)
                {
                  aij_bc_value = 0.0;
                  for (size_t qp : fe_srf_data.QuadraturePointIndices())
                    aij_bc_value += kappa * source_function_->Evaluate(fe_srf_data.QPointXYZ(qp)) *
                                    fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(jm, qp) *
                                    fe_srf_data.JxW(qp);
                }

                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for fi

            // Assemble gradient terms
            // For the following comments we use the notation:
            // Dk = 0.5* n dot nabla bk

            // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
            for (size_t i = 0; i < num_nodes; ++i)
            {
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t j = 0; j < num_nodes; ++j)
              {
                Vector3 vec_aij;
                for (size_t qp : fe_srf_data.QuadraturePointIndices())
                  vec_aij += fe_srf_data.ShapeValue(j, qp) * fe_srf_data.ShapeGrad(i, qp) *
                               fe_srf_data.JxW(qp) +
                             fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeGrad(j, qp) *
                               fe_srf_data.JxW(qp);
                const double aij = -Dg * n_f.Dot(vec_aij);

                double aij_bc_value = aij * bc_value;

                if (ref_solution_function_)
                {
                  Vector3 vec_aij_mms;
                  for (size_t qp : fe_srf_data.QuadraturePointIndices())
                    vec_aij_mms += ref_solution_function_->Evaluate(fe_srf_data.QPointXYZ(qp)) *
                                   (fe_srf_data.ShapeValue(j, qp) * fe_srf_data.ShapeGrad(i, qp) *
                                      fe_srf_data.JxW(qp) +
                                    fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeGrad(j, qp) *
                                      fe_srf_data.JxW(qp));
                  aij_bc_value = -Dg * n_f.Dot(vec_aij_mms);
                }

                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for i
          }     // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12)
              continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(fval) >= 1.0e-12)
              {
                double rhs_val = 0.0;
                for (size_t qp : fe_srf_data.QuadraturePointIndices())
                  rhs_val += fe_srf_data.ShapeValue(i, qp) * fe_srf_data.JxW(qp);
                rhs_val *= (fval / bval);

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

  KSPSetOperators(ksp_, A_, A_);

  if (options.verbose)
    log.Log() << program_timer.GetTimeString() << " Assembly completed";

  PC pc;
  KSPGetPC(ksp_, &pc);
  PCSetUp(pc);

  KSPSetUp(ksp_);
}

void
DiffusionMIPSolver::AssembleAand_b(const std::vector<double>& q_vector)
{
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

    for (size_t g = 0; g < num_groups; ++g)
    {
      // Get coefficient and nodal src
      const double Dg = xs.Dg[g];
      const double sigr_g = xs.sigR[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; ++j)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      // Assemble continuous terms
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0;
        for (size_t j = 0; j < num_nodes; ++j)
        {
          const int64_t jmap = sdm_.MapDOF(cell, j, uk_man_, 0, g);

          const double entry_aij =
            Dg * intV_gradshapeI_gradshapeJ[i][j] + sigr_g * intV_shapeI_shapeJ[i][j];

          entry_rhs_i += intV_shapeI_shapeJ[i][j] * qg[j];

          MatSetValue(A_, imap, jmap, entry_aij, ADD_VALUES);
        } // for j

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      // Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const auto& n_f = face.normal;
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& intS_shapeI_shapeJ = unit_cell_matrices.intS_shapeI_shapeJ[f];
        const auto& intS_shapeI_gradshapeJ = unit_cell_matrices.intS_shapeI_gradshapeJ[f];
        const auto& intS_shapeI = unit_cell_matrices.intS_shapeI[f];

        const double hm = HPerpendicular(cell, f);

        if (face.has_neighbor)
        {
          const auto& adj_cell = grid_.cells[face.neighbor_id];
          const auto& adj_cell_mapping = sdm_.GetCellMapping(adj_cell);
          const auto ac_nodes = adj_cell_mapping.GetNodeLocations();
          const size_t acf = MeshContinuum::MapCellFace(cell, adj_cell, f);
          const double hp = HPerpendicular(adj_cell, acf);

          const auto& adj_xs = mat_id_2_xs_map_.at(adj_cell.material_id);
          const double adj_Dg = adj_xs.Dg[g];

          // Compute kappa
          double kappa = 1.0;
          if (cell.Type() == CellType::SLAB)
            kappa = fmax(options.penalty_factor * (adj_Dg / hp + Dg / hm) * 0.5, 0.25);
          if (cell.Type() == CellType::POLYGON)
            kappa = fmax(options.penalty_factor * (adj_Dg / hp + Dg / hm) * 0.5, 0.25);
          if (cell.Type() == CellType::POLYHEDRON)
            kappa = fmax(options.penalty_factor * 2.0 * (adj_Dg / hp + Dg / hm) * 0.5, 0.25);

          // Assembly penalty terms
          for (size_t fi = 0; fi < num_face_nodes; ++fi)
          {
            const int i = cell_mapping.MapFaceNode(f, fi);
            const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

            for (size_t fj = 0; fj < num_face_nodes; ++fj)
            {
              const int jm = cell_mapping.MapFaceNode(f, fj); // j-minus
              const int jp =
                MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fj); // j-plus
              const int64_t jmmap = sdm_.MapDOF(cell, jm, uk_man_, 0, g);
              const int64_t jpmap = sdm_.MapDOF(adj_cell, jp, uk_man_, 0, g);

              const double aij = kappa * intS_shapeI_shapeJ[i][jm];

              MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
              MatSetValue(A_, imap, jpmap, -aij, ADD_VALUES);
            } // for fj
          }   // for fi

          // Assemble gradient terms
          // For the following comments we use the notation:
          // Dk = 0.5* n dot nabla bk

          // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
          for (int i = 0; i < num_nodes; ++i)
          {
            const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

            for (int fj = 0; fj < num_face_nodes; ++fj)
            {
              const int jm = cell_mapping.MapFaceNode(f, fj); // j-minus
              const int jp =
                MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fj); // j-plus
              const int64_t jmmap = sdm_.MapDOF(cell, jm, uk_man_, 0, g);
              const int64_t jpmap = sdm_.MapDOF(adj_cell, jp, uk_man_, 0, g);

              const double aij = -0.5 * Dg * n_f.Dot(intS_shapeI_gradshapeJ[jm][i]);

              MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
              MatSetValue(A_, imap, jpmap, -aij, ADD_VALUES);
            } // for fj
          }   // for i

          // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
          for (int fi = 0; fi < num_face_nodes; ++fi)
          {
            const int im = cell_mapping.MapFaceNode(f, fi); // i-minus
            const int ip =
              MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fi); // i-plus
            const int64_t immap = sdm_.MapDOF(cell, im, uk_man_, 0, g);
            const int64_t ipmap = sdm_.MapDOF(adj_cell, ip, uk_man_, 0, g);

            for (int j = 0; j < num_nodes; ++j)
            {
              const int64_t jmap = sdm_.MapDOF(cell, j, uk_man_, 0, g);

              const double aij = -0.5 * Dg * n_f.Dot(intS_shapeI_gradshapeJ[im][j]);

              MatSetValue(A_, immap, jmap, aij, ADD_VALUES);
              MatSetValue(A_, ipmap, jmap, -aij, ADD_VALUES);
            } // for j
          }   // for fi

        } // internal face
        else
        {
          BoundaryCondition bc;
          if (bcs_.count(face.neighbor_id) > 0)
            bc = bcs_.at(face.neighbor_id);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            // Compute kappa
            double kappa = 1.0;
            if (cell.Type() == CellType::SLAB)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYGON)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYHEDRON)
              kappa = fmax(options.penalty_factor * 2.0 * Dg / hm, 0.25);

            // Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const int jm = cell_mapping.MapFaceNode(f, fj);
                const int64_t jmmap = sdm_.MapDOF(cell, jm, uk_man_, 0, g);

                const double aij = kappa * intS_shapeI_shapeJ[i][jm];
                const double aij_bc_value = aij * bc_value;

                MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for fi

            // Assemble gradient terms
            // For the following comments we use the notation:
            // Dk = n dot nabla bk

            // D* n dot (b_j^+ - b_j^-)*nabla b_i^-
            for (size_t i = 0; i < num_nodes; ++i)
            {
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t j = 0; j < num_nodes; ++j)
              {
                const int64_t jmap = sdm_.MapDOF(cell, j, uk_man_, 0, g);

                const double aij =
                  -Dg * n_f.Dot(intS_shapeI_gradshapeJ[j][i] + intS_shapeI_gradshapeJ[i][j]);
                const double aij_bc_value = aij * bc_value;

                MatSetValue(A_, imap, jmap, aij, ADD_VALUES);
                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for i
          }     // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double aval = bc.values[0];
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12)
              continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(aval) >= 1.0e-12)
              {
                for (size_t fj = 0; fj < num_face_nodes; ++fj)
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
DiffusionMIPSolver::Assemble_b(const std::vector<double>& q_vector)
{
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

    const auto& intV_shapeI_shapeJ = unit_cell_matrices.intV_shapeI_shapeJ;

    const auto& xs = mat_id_2_xs_map_.at(cell.material_id);

    for (size_t g = 0; g < num_groups; ++g)
    {
      // Get coefficient and nodal src
      const double Dg = xs.Dg[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; ++j)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      // Assemble continuous terms
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0;
        for (size_t j = 0; j < num_nodes; ++j)
          entry_rhs_i += intV_shapeI_shapeJ[i][j] * qg[j];

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      // Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const auto& n_f = face.normal;
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& intS_shapeI_shapeJ = unit_cell_matrices.intS_shapeI_shapeJ[f];
        const auto& intS_shapeI_gradshapeJ = unit_cell_matrices.intS_shapeI_gradshapeJ[f];
        const auto& intS_shapeI = unit_cell_matrices.intS_shapeI[f];

        const double hm = HPerpendicular(cell, f);

        if (not face.has_neighbor and not suppress_bcs_)
        {
          BoundaryCondition bc;
          if (bcs_.count(face.neighbor_id) > 0)
            bc = bcs_.at(face.neighbor_id);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            // Compute kappa
            double kappa = 1.0;
            if (cell.Type() == CellType::SLAB)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYGON)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYHEDRON)
              kappa = fmax(options.penalty_factor * 2.0 * Dg / hm, 0.25);

            // Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const int jm = cell_mapping.MapFaceNode(f, fj);

                const double aij = kappa * intS_shapeI_shapeJ[i][jm];
                const double aij_bc_value = aij * bc_value;

                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for fi

            // Assemble gradient terms
            // For the following comments we use the notation:
            // Dk = n dot nabla bk

            // D* n dot (b_j^+ - b_j^-)*nabla b_i^-
            for (size_t i = 0; i < num_nodes; ++i)
            {
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t j = 0; j < num_nodes; ++j)
              {
                const double aij =
                  -Dg * n_f.Dot(intS_shapeI_gradshapeJ[j][i] + intS_shapeI_gradshapeJ[i][j]);
                const double aij_bc_value = aij * bc_value;

                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for i
          }     // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12)
              continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; ++fi)
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
DiffusionMIPSolver::Assemble_b(Vec petsc_q_vector)
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

    const auto& xs = mat_id_2_xs_map_.at(cell.material_id);

    for (size_t g = 0; g < num_groups; ++g)
    {
      // Get coefficient and nodal src
      const double Dg = xs.Dg[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; ++j)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      // Assemble continuous terms
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0;
        for (size_t j = 0; j < num_nodes; ++j)
          entry_rhs_i += intV_shapeI_shapeJ[i][j] * qg[j];

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      // Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        const auto& n_f = face.normal;
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);

        const auto& intS_shapeI_shapeJ = unit_cell_matrices.intS_shapeI_shapeJ[f];
        const auto& intS_shapeI_gradshapeJ = unit_cell_matrices.intS_shapeI_gradshapeJ[f];
        const auto& intS_shapeI = unit_cell_matrices.intS_shapeI[f];

        const double hm = HPerpendicular(cell, f);

        if (not face.has_neighbor and not suppress_bcs_)
        {
          BoundaryCondition bc;
          if (bcs_.count(face.neighbor_id) > 0)
            bc = bcs_.at(face.neighbor_id);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            // Compute kappa
            double kappa = 1.0;
            if (cell.Type() == CellType::SLAB)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYGON)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == CellType::POLYHEDRON)
              kappa = fmax(options.penalty_factor * 2.0 * Dg / hm, 0.25);

            // Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const int jm = cell_mapping.MapFaceNode(f, fj);

                const double aij = kappa * intS_shapeI_shapeJ[i][jm];
                const double aij_bc_value = aij * bc_value;

                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for fi

            // Assemble gradient terms
            // For the following comments we use the notation:
            // Dk = n dot nabla bk

            // D* n dot (b_j^+ - b_j^-)*nabla b_i^-
            for (size_t i = 0; i < num_nodes; ++i)
            {
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t j = 0; j < num_nodes; ++j)
              {
                const double aij =
                  -Dg * n_f.Dot(intS_shapeI_gradshapeJ[j][i] + intS_shapeI_gradshapeJ[i][j]);
                const double aij_bc_value = aij * bc_value;

                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for i
          }     // Dirichlet BC
          else if (bc.type == BCType::ROBIN)
          {
            const double bval = bc.values[1];
            const double fval = bc.values[2];

            if (std::fabs(bval) < 1.0e-12)
              continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; ++fi)
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

double
DiffusionMIPSolver::HPerpendicular(const Cell& cell, unsigned int f)
{
  const auto& cell_mapping = sdm_.GetCellMapping(cell);
  double hp;

  const size_t num_faces = cell.faces.size();
  const size_t num_vertices = cell.vertex_ids.size();

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
  if (cell.Type() == CellType::SLAB)
    hp = volume / 2.0;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell.Type() == CellType::POLYGON)
  {
    if (num_faces == 3)
      hp = 2.0 * volume / face_area;
    else if (num_faces == 4)
      hp = volume / face_area;
    else // Nv > 4
    {
      const double surface_area = ComputeSurfaceArea();

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
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell.Type() == CellType::POLYHEDRON)
  {
    const double surface_area = ComputeSurfaceArea();

    if (num_faces == 4) // Tet
      hp = 3 * volume / surface_area;
    else if (num_faces == 6 and num_vertices == 8) // Hex
      hp = volume / surface_area;
    else // Polyhedron
      hp = 6 * volume / surface_area;
  } // Polyhedron
  else
    throw std::logic_error("acceleration::DiffusionMIPSolver::HPerpendicular: "
                           "Unsupported cell type in call to HPerpendicular");

  return hp;
}

int
DiffusionMIPSolver::MapFaceNodeDisc(const Cell& cur_cell,
                                    const Cell& adj_cell,
                                    const std::vector<Vector3>& cc_node_locs,
                                    const std::vector<Vector3>& ac_node_locs,
                                    size_t ccf,
                                    size_t acf,
                                    size_t ccfi,
                                    double epsilon)
{
  const auto& cur_cell_mapping = sdm_.GetCellMapping(cur_cell);
  const auto& adj_cell_mapping = sdm_.GetCellMapping(adj_cell);

  const int i = cur_cell_mapping.MapFaceNode(ccf, ccfi);
  const auto& node_i_loc = cc_node_locs[i];

  const size_t adj_face_num_nodes = adj_cell_mapping.NumFaceNodes(acf);

  for (size_t fj = 0; fj < adj_face_num_nodes; ++fj)
  {
    const int j = adj_cell_mapping.MapFaceNode(acf, fj);
    if ((node_i_loc - ac_node_locs[j]).NormSquare() < epsilon)
      return j;
  }

  throw std::logic_error("acceleration::DiffusionMIPSolver::MapFaceNodeDisc: Mapping failure.");
}

} // namespace opensn
