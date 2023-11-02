#include "chi_ffinter_line.h"

#include "mesh/Cell/cell.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Initializes the data structures necessary for interpolation. This is
 * independent of the physics and hence is a routine on its own.*/
void chi_mesh::FieldFunctionInterpolationLine::
  Initialize()
{
  Chi::log.Log0Verbose1() << "Initializing line interpolator.";
  //================================================== Check for empty FF-list
  if (field_functions_.empty())
    throw std::logic_error("Unassigned field function in line "
                           "field function interpolator.");

  //================================================== Create points;
  const chi_mesh::Vector3 vif = pf_ - pi_;
  delta_d_ = vif.Norm() / (number_of_points_ - 1);

  const auto omega = vif.Normalized();

  interpolation_points_.push_back(pi_);
  for (int k=1; k<(number_of_points_); k++)
    interpolation_points_.push_back(pi_ + omega * delta_d_ * k);

  //====================================================== Loop over contexts
  const size_t num_ff = field_functions_.size();
  for (size_t ff=0; ff<num_ff; ff++)
  {
    ff_contexts_.emplace_back();
    auto& ff_context = ff_contexts_.back();

    ff_context.ref_ff = field_functions_[ff];
    const auto& sdm = ff_context.ref_ff->GetSpatialDiscretization();
    const auto& grid = sdm.Grid();

    ff_context.interpolation_points_ass_cell.assign(number_of_points_, 0);
    ff_context.interpolation_points_has_ass_cell.assign(number_of_points_, false);

    //================================================== Find a home for each
    //                                                   point
    for (const auto& cell : grid.local_cells)
    {
      for (int p=0; p < number_of_points_; p++)
      {
        const auto& point = interpolation_points_[p];
        if (grid.CheckPointInsideCell(cell, point))
        {
          ff_context.interpolation_points_ass_cell[p] = cell.local_id_;
          ff_context.interpolation_points_has_ass_cell[p] = true;
        }
      }//for point p
    }//for cell
  }//for ff

  Chi::log.Log0Verbose1() << "Finished initializing interpolator.";
}