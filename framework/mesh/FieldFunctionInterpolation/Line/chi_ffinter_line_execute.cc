#include "chi_ffinter_line.h"

#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Executes the interpolation.*/
void chi_mesh::FieldFunctionInterpolationLine::Execute()
{
  Chi::log.Log0Verbose1() << "Executing line interpolator.";
  for (int ff=0; ff < field_functions_.size(); ff++)
  {
          auto& ff_ctx = ff_contexts_[ff];
    const auto& ref_ff = *ff_ctx.ref_ff;
    const auto& sdm    = ref_ff.GetSpatialDiscretization();
    const auto& grid   = sdm.Grid();

    const auto& uk_man = ref_ff.GetUnknownManager();
    const auto uid = 0;
    const auto cid = ref_component_;

    const auto field_data = ref_ff.GetGhostedFieldVector();

    ff_ctx.interpolation_points_values.assign(number_of_points_, 0.0);
    for (int p=0; p < number_of_points_; ++p)
    {
      if (not ff_ctx.interpolation_points_has_ass_cell[p]) continue;

      const auto cell_local_index = ff_ctx.interpolation_points_ass_cell[p];
      const auto& cell = grid.local_cells[cell_local_index];
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      std::vector<double> shape_function_vals(num_nodes, 0.0);
      cell_mapping.ShapeValues(interpolation_points_[p], shape_function_vals);

      double point_value = 0.0;
      for (size_t i=0; i<num_nodes; ++i)
      {
        const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);

        point_value += shape_function_vals[i]*field_data[imap];
      }//for node i
      ff_ctx.interpolation_points_values[p] = point_value;
    }//for p
  }//for ff

}