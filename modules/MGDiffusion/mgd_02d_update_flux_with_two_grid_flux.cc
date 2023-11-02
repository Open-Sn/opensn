#include "mg_diffusion_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

//============================================= assemble matrix A
void
mg_diffusion::Solver::Update_Flux_With_TwoGrid(const int64_t verbose)
{
  if (verbose > 2) Chi::log.Log() << "\nUpdating Thermal fluxes from two-grid";

  const auto& grid = *grid_ptr_;
  const auto& sdm = *sdm_ptr_;

  // contains two_grid flux, stored in last num_groups entry
  const double* xlocal_tg;
  VecGetArrayRead(x_[num_groups_], &xlocal_tg);

  int counter = 0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& xstg = map_mat_id_2_tginfo.at(cell.material_id_);

    for (unsigned int g = last_fast_group_; g < num_groups_; ++g)
    {
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imap = sdm.MapDOFLocal(cell, i);
        const int64_t imap_globl = sdm.MapDOF(cell, i);
        const double aux = xlocal_tg[imap] * VF_[counter][i] * xstg.spectrum[g];
        VecSetValue(x_[g], imap_globl, aux, ADD_VALUES);
      } // i
    }   // g
    counter++;
  } // for cell

  // finalize
  for (unsigned int g = last_fast_group_; g < num_groups_; ++g)
  {
    VecAssemblyBegin(x_[g]);
    VecAssemblyEnd(x_[g]);
  }
  // release two-grid flux
  VecRestoreArrayRead(x_[num_groups_], &xlocal_tg);
}