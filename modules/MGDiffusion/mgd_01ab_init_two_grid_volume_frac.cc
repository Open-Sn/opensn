#include "mg_diffusion_solver.h"
#include "utils/chi_timer.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//============================================= assemble matrix A
void mg_diffusion::Solver::Compute_TwoGrid_VolumeFractions()
{
  const auto& grid = *grid_ptr_;
  const auto& sdm  = *sdm_ptr_;

  const size_t ncells = grid.local_cells.size();
  VF_.resize(ncells);

  int counter = 0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto  qp_data      = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes   = cell_mapping.NumNodes();

    VF_[counter].resize(num_nodes, 0.0);

    for (size_t i=0; i < num_nodes; ++i)
    {
      double vol_frac_shape_i = 0.0;
      for (size_t qp : qp_data.QuadraturePointIndices())
        vol_frac_shape_i +=  qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      vol_frac_shape_i /= cell_mapping.CellVolume();
      VF_[counter][i] = vol_frac_shape_i;
    }//for i

    counter++;
  }//for cell

}
