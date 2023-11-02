#include "fieldfunction_gridbased.h"

// #########################################################
/**Makes a ghosted version of the field vector.*/
std::vector<double>
chi_physics::FieldFunctionGridBased::GetGhostedFieldVector() const
{
  return ghosted_field_vector_->LocalSTLData();
}