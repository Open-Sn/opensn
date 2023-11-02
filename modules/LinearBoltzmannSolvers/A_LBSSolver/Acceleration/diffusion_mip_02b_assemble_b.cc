#include "diffusion_mip.h"
#include "acceleration.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

#include "physics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"

#include "chi_runtime.h" //TODO:Remove
#include "chi_log.h"     //TODO:Remove
#include "utils/chi_timer.h"
#include "console/chi_console.h"

#define DefaultBCDirichlet                                                                         \
  BoundaryCondition                                                                                \
  {                                                                                                \
    BCType::DIRICHLET,                                                                             \
    {                                                                                              \
      0, 0, 0                                                                                      \
    }                                                                                              \
  }

//###################################################################
/**Assembles just the RHS using quadrature points. These
 * routines exist for implementing MMS.*/
void
lbs::acceleration::DiffusionMIPSolver::Assemble_b_wQpoints(const std::vector<double>& q_vector)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::"
                            "AssembleAand_b_wQpoints";
  if (A_ == nullptr or rhs_ == nullptr or ksp_ == nullptr)
    throw std::logic_error(fname + ": Some or all PETSc elements are null. "
                                   "Check that Initialize has been called.");
  if (options.verbose) Chi::log.Log() << Chi::program_timer.GetTimeString() << " Starting assembly";

  lua_State* L = Chi::console.GetConsoleState();
  const auto& source_function = options.source_lua_function;
  const auto& solution_function = options.ref_solution_lua_function;

  VecSet(rhs_, 0.0);

  for (const auto& cell : grid_.local_cells)
  {
    const size_t num_faces = cell.faces_.size();
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_groups = uk_man_.unknowns_.front().num_components_;

    const auto& xs = mat_id_2_xs_map_.at(cell.material_id_);

    //=========================================== For component/group
    for (size_t g = 0; g < num_groups; ++g)
    {
      //==================================== Get coefficient and nodal src
      const double Dg = xs.Dg[g];

      std::vector<double> qg(num_nodes, 0.0);
      for (size_t j = 0; j < num_nodes; j++)
        qg[j] = q_vector[sdm_.MapDOFLocal(cell, j, uk_man_, 0, g)];

      //==================================== Assemble continuous terms
      for (size_t i = 0; i < num_nodes; i++)
      {
        const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);
        double entry_rhs_i = 0.0; // entry may accumulate over j
        if (source_function.empty())
          for (size_t j = 0; j < num_nodes; j++)
          {
            for (size_t qp : qp_data.QuadraturePointIndices())
            {
              entry_rhs_i +=
                qg[j] * qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp) * qp_data.JxW(qp);
            } // for qp
          }   // for j
        else
        {
          for (size_t qp : qp_data.QuadraturePointIndices())
            entry_rhs_i += CallLuaXYZFunction(L, source_function, qp_data.QPointXYZ(qp)) *
                           qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
        }

        VecSetValue(rhs_, imap, entry_rhs_i, ADD_VALUES);
      } // for i

      //==================================== Assemble face terms
      for (size_t f = 0; f < num_faces; ++f)
      {
        const auto& face = cell.faces_[f];
        const auto& n_f = face.normal_;
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        const auto fqp_data = cell_mapping.MakeSurfaceQuadraturePointData(f);

        const double hm = HPerpendicular(cell, f);

        if (not face.has_neighbor_)
        {
          auto bc = DefaultBCDirichlet;
          if (bcs_.count(face.neighbor_id_) > 0) bc = bcs_.at(face.neighbor_id_);

          if (bc.type == BCType::DIRICHLET)
          {
            const double bc_value = bc.values[0];

            //========================= Compute kappa
            double kappa = 1.0;
            if (cell.Type() == chi_mesh::CellType::SLAB)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == chi_mesh::CellType::POLYGON)
              kappa = fmax(options.penalty_factor * Dg / hm, 0.25);
            if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
              kappa = fmax(options.penalty_factor * 2.0 * Dg / hm, 0.25);

            //========================= Assembly penalty terms
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const int jm = cell_mapping.MapFaceNode(f, fj);

                double aij = 0.0;
                for (size_t qp : fqp_data.QuadraturePointIndices())
                  aij += kappa * fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) *
                         fqp_data.JxW(qp);
                double aij_bc_value = aij * bc_value;

                if (not solution_function.empty())
                {
                  aij_bc_value = 0.0;
                  for (size_t qp : fqp_data.QuadraturePointIndices())
                    aij_bc_value +=
                      kappa * CallLuaXYZFunction(L, solution_function, fqp_data.QPointXYZ(qp)) *
                      fqp_data.ShapeValue(i, qp) * fqp_data.ShapeValue(jm, qp) * fqp_data.JxW(qp);
                }

                VecSetValue(rhs_, imap, aij_bc_value, ADD_VALUES);
              } // for fj
            }   // for fi

            //========================= Assemble gradient terms
            // For the following comments we use the notation:
            // Dk = 0.5* n dot nabla bk

            // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
            for (size_t i = 0; i < num_nodes; i++)
            {
              const int64_t imap = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              for (size_t j = 0; j < num_nodes; j++)
              {
                chi_mesh::Vector3 vec_aij;
                for (size_t qp : fqp_data.QuadraturePointIndices())
                  vec_aij +=
                    fqp_data.ShapeValue(j, qp) * fqp_data.ShapeGrad(i, qp) * fqp_data.JxW(qp) +
                    fqp_data.ShapeValue(i, qp) * fqp_data.ShapeGrad(j, qp) * fqp_data.JxW(qp);
                const double aij = -Dg * n_f.Dot(vec_aij);

                double aij_bc_value = aij * bc_value;

                if (not solution_function.empty())
                {
                  chi_mesh::Vector3 vec_aij_mms;
                  for (size_t qp : fqp_data.QuadraturePointIndices())
                    vec_aij_mms +=
                      CallLuaXYZFunction(L, solution_function, fqp_data.QPointXYZ(qp)) *
                      (fqp_data.ShapeValue(j, qp) * fqp_data.ShapeGrad(i, qp) * fqp_data.JxW(qp) +
                       fqp_data.ShapeValue(i, qp) * fqp_data.ShapeGrad(j, qp) * fqp_data.JxW(qp));
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

            if (std::fabs(bval) < 1.0e-12) continue; // a and f assumed zero

            for (size_t fi = 0; fi < num_face_nodes; fi++)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const int64_t ir = sdm_.MapDOF(cell, i, uk_man_, 0, g);

              if (std::fabs(fval) >= 1.0e-12)
              {
                double rhs_val = 0.0;
                for (size_t qp : fqp_data.QuadraturePointIndices())
                  rhs_val += fqp_data.ShapeValue(i, qp) * fqp_data.JxW(qp);
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
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " Assembly completed";

  PC pc;
  KSPGetPC(ksp_, &pc);
  PCSetUp(pc);

  KSPSetUp(ksp_);
}