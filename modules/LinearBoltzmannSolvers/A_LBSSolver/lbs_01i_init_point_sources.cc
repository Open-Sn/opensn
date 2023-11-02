#include "lbs_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

/**Intializes all point sources.*/
void lbs::LBSSolver::InitializePointSources()
{
  const std::string fname = "InitializePointSources";

  //============================================= Loop over point sources
  for (auto& point_source : point_sources_)
  {
    if (point_source.Strength().size() != num_groups_)
      throw std::logic_error(
        fname + ": Point source multigroup strength vector "
                "is not compatible with the number of "
                "groups in the simulation. Expected " +
        std::to_string(num_groups_) + " found " +
        std::to_string(point_source.Strength().size()));

    const auto& p = point_source.Location();
    double v_total = 0.0; //Total volume of all cells sharing
                          // this source
    std::vector<PointSource::ContainingCellInfo> temp_list;
    for (const auto& cell : grid_ptr_->local_cells)
    {
      if (grid_ptr_->CheckPointInsideCell(cell, p))
      {
        const auto& cell_view = discretization_->GetCellMapping(cell);
        const auto& cell_matrices = unit_cell_matrices_[cell.local_id_];
        const auto& M = cell_matrices.M_matrix;
        const auto& I = cell_matrices.Vi_vectors;

        std::vector<double> shape_values;
        cell_view.ShapeValues(point_source.Location(),
                              shape_values/**ByRef*/);

        const auto M_inv = chi_math::Inverse(M);

        const auto q_p_weights = chi_math::MatMul(M_inv, shape_values);

        double v_cell = 0.0;
        for (double val : I) v_cell += val;
        v_total += v_cell;

        temp_list.push_back(
          PointSource::ContainingCellInfo{v_cell,
                                          cell.local_id_,
                                          shape_values,
                                          q_p_weights});
      }//if inside
    }//for local cell

    auto ghost_global_ids = grid_ptr_->cells.GetGhostGlobalIDs();
    for (uint64_t ghost_global_id : ghost_global_ids)
    {
      const auto& neighbor_cell = grid_ptr_->cells[ghost_global_id];
      if (grid_ptr_->CheckPointInsideCell(neighbor_cell, p))
      {
        const auto& cell_matrices =
          unit_ghost_cell_matrices_[neighbor_cell.global_id_];
        for (double val : cell_matrices.Vi_vectors)
          v_total += val;
      }//if point inside
    }//for ghost cell

    point_source.ClearInitializedInfo();
    for (const auto& info : temp_list)
    {
      point_source.AddContainingCellInfo(info.volume_weight/v_total,
                                         info.cell_local_id,
                                         info.shape_values,
                                         info.node_weights);
      const auto& cell = grid_ptr_->local_cells[info.cell_local_id];
      //Output message
      {
        std::stringstream output;
        output << "Point source at " << p.PrintStr() << " assigned to cell "
               << cell.global_id_ << " with shape values ";
        for (double val : info.shape_values) output << val << " ";
        output << "volume_weight=" << info.volume_weight/v_total;

        Chi::log.LogAll() << output.str();
      }
    }//for info in temp list
  }//for point_source
}