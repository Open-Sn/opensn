#include "chi_ffinter_volume.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Initializes the volume field function interpolation.*/
void chi_mesh::FieldFunctionInterpolationVolume::Initialize()
{
  Chi::log.Log0Verbose1() << "Initializing volume interpolator.";
  //================================================== Check grid available
  if (field_functions_.empty())
    throw std::logic_error("Unassigned field function in volume field "
                           "function interpolator.");

  if (logical_volume_ == nullptr)
    throw std::logic_error("Unassigned logical volume in volume field function"
                           "interpolator.");

  const auto& grid =
    field_functions_.front()->GetSpatialDiscretization().Grid();

  //================================================== Find cells inside volume
  for (const auto& cell : grid.local_cells)
    if (logical_volume_->Inside(cell.centroid_))
      cell_local_ids_inside_logvol_.push_back(cell.local_id_);
}