#include "fieldfunction_gridbased.h"

#include <petsc.h>

#include "chi_log.h"

// ###################################################################
/**Updates the field data with a STL vector.*/
void chi_physics::FieldFunctionGridBased::UpdateFieldVector(
  const std::vector<double>& field_vector)
{
  ChiInvalidArgumentIf(field_vector.size() < ghosted_field_vector_->LocalSize(),
                       "Attempted update with a vector of insufficient size.");

  ghosted_field_vector_->Set(field_vector);

  ghosted_field_vector_->CommunicateGhostEntries();
}

// ###################################################################
/**Updates the field data with a PETSc vector.*/
void chi_physics::FieldFunctionGridBased::UpdateFieldVector(
  const Vec& field_vector)
{
  ghosted_field_vector_->CopyLocalValues(field_vector);

  ghosted_field_vector_->CommunicateGhostEntries();
}