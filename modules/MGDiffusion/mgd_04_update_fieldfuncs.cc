#include "mg_diffusion_solver.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"

//###################################################################
/**Updates the field functions with the lates data.*/
void mg_diffusion::Solver::UpdateFieldFunctions()
{
  const auto& OneDOFPerNode = sdm_ptr_->UNITARY_UNKNOWN_MANAGER;
  for (int g=0; g < num_groups_; ++g)
  {
    std::vector<double> data_vector;
    sdm_ptr_->LocalizePETScVector(x_[g], data_vector, OneDOFPerNode);

    auto& ff = field_functions_.at(g);
    ff->UpdateFieldVector(data_vector);
  }//for g
}