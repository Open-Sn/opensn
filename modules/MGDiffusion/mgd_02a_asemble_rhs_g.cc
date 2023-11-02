#include "mg_diffusion_solver.h"
#include "chi_runtime.h"
#include "chi_log.h"
#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//========================================================== Solve 1g problem
void mg_diffusion::Solver::Assemble_RHS(const unsigned int g,
                                        const int64_t verbose)
{
  if (verbose > 2)
    Chi::log.Log() << "\nAssemblying RHS for group " + std::to_string(g);

  // copy the external source vector for group g into b
  VecSet(b_, 0.0);
  VecCopy(bext_[g], b_);

  const auto& sdm  = *mg_diffusion::Solver::sdm_ptr_;
  // compute inscattering term
  for (const auto& cell :  mg_diffusion::Solver::grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto  qp_data      = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes   = cell_mapping.NumNodes();

    const auto& xs = matid_to_xs_map.at(cell.material_id_);
    const auto& S = xs->TransferMatrix(0);

    for (const auto& [row_g, gprime, sigma_sm] : S.Row(g))
    {
      if (gprime != g) // g and row_g are the same, maybe different int types
      {
        const double* xlocal;
        VecGetArrayRead(x_[gprime], &xlocal);

        for (size_t i=0; i<num_nodes; ++i)
        {
          const int64_t imap = sdm.MapDOF(cell,i);
          double inscatter_g = 0.0;

          for (size_t j = 0; j < num_nodes; ++j)
          {
            const int64_t jmap = sdm.MapDOFLocal(cell, j);

           // get flux at node j
            const double flxj_gp = xlocal[jmap];
            for (size_t qp: qp_data.QuadraturePointIndices())
              inscatter_g += sigma_sm * flxj_gp *
                         qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp) *
                         qp_data.JxW(qp);
          }//for j
          // add inscattering value to vector
          VecSetValue(b_, imap, inscatter_g, ADD_VALUES);
        }//for i
        VecRestoreArrayRead(x_[gprime], &xlocal);
      }//if gp!=g
    }// for gprime
  }//for cell

  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);

}

//  cout << "b=bext["<<g<<"]+inscattering\n";
//  VecView(b, PETSC_VIEWER_STDERR_WORLD);
