#include "petsc_utils.h"

#include "math/ParallelVector/ParallelVector.h"

#include "chi_log.h"

// ###################################################################
/**Copies a PETSc vector to a STL vector. Only the local portion is
 * copied.*/
void chi_math::PETScUtils::CopyVecToSTLvector(Vec x,
                                              std::vector<double>& data,
                                              size_t N,
                                              bool resize_STL /*=true*/)
{
  if (resize_STL)
  {
    data.clear();
    data.assign(N, 0.0);
  }
  else
    ChiLogicalErrorIf(data.size() < N,
                      "data.size() < N, " + std::to_string(data.size()) +
                        " < " + std::to_string(N));

  const double* x_ref;
  VecGetArrayRead(x, &x_ref);

  std::copy(x_ref, x_ref + N, data.begin());

  VecRestoreArrayRead(x, &x_ref);
}

// ###################################################################
/**Copies a PETSc vector to a STL vector. Only the local portion is
 * copied.*/
void chi_math::PETScUtils::CopyVecToSTLvectorWithGhosts(
  Vec x, std::vector<double>& data, size_t N, bool resize_STL /*=true*/)
{
  if (resize_STL)
  {
    data.clear();
    data.assign(N, 0.0);
  }
  else
    ChiLogicalErrorIf(data.size() != N,
                      "data.size() != N, " + std::to_string(data.size()) +
                        " < " + std::to_string(N));

  auto info = GetGhostVectorLocalViewRead(x);
  const double* x_ref = info.x_localized_raw;

  std::copy(x_ref, x_ref + N, data.begin());

  RestoreGhostVectorLocalViewRead(x, info);
}

void chi_math::PETScUtils::CopySTLvectorToVec(const std::vector<double>& data,
                                              Vec x,
                                              size_t N)
{
  double* x_ref;
  VecGetArray(x, &x_ref);

  std::copy(data.begin(), data.end(), x_ref);

  VecRestoreArray(x, &x_ref);
}

void chi_math::PETScUtils::CopyParallelVectorToVec(const ParallelVector& y,
                                                   Vec x)
{
  const double* y_data = y.Data();
  double* x_data;
  VecGetArray(x, &x_data);
  std::copy(y_data, y_data + y.LocalSize(), x_data);
  VecRestoreArray(x, &x_data);
}

// ###################################################################
/**Copies global values from a PETSc vector to a STL vector.*/
void chi_math::PETScUtils::CopyGlobalVecToSTLvector(
  Vec x, const std::vector<int64_t>& global_indices, std::vector<double>& data)
{
  //=================================== Populating local indices
  size_t N = global_indices.size();
  std::vector<int64_t> local_indices(N, 0);
  size_t counter = 0;
  for (auto val : global_indices)
  {
    local_indices[counter] = counter;
    ++counter;
  }

  //=================================== Creating PETSc vector
  Vec local_vec;
  VecCreateSeq(PETSC_COMM_SELF, global_indices.size() + 1, &local_vec);
  VecSet(local_vec, 0.0);

  //=================================== Create and transfer index sets
  IS global_set;
  IS local_set;
  ISCreateGeneral(
    PETSC_COMM_SELF, N, global_indices.data(), PETSC_COPY_VALUES, &global_set);
  ISCreateGeneral(
    PETSC_COMM_SELF, N, local_indices.data(), PETSC_COPY_VALUES, &local_set);
  VecScatter scat;
  VecScatterCreate(x, global_set, local_vec, local_set, &scat);
  VecScatterBegin(scat, x, local_vec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scat, x, local_vec, INSERT_VALUES, SCATTER_FORWARD);

  //=================================== Copy to STL
  data.clear();
  data.resize(N, 0.0);
  const double* x_ref;
  VecGetArrayRead(local_vec, &x_ref);

  std::copy(x_ref, x_ref + N, data.begin());

  VecRestoreArrayRead(x, &x_ref);

  //=================================== Cleanup
  ISDestroy(&global_set);
  ISDestroy(&local_set);

  VecDestroy(&local_vec);
}

// ###################################################################
/**Communicates ghost entries of a ghost vector. This operation
 * is suitable when only a single vector is communicated. When
 * more than vector is communicated it would be more efficient
 * to "Begin" all the vectors followed by and "End" of each
 * vector.*/
void chi_math::PETScUtils::CommunicateGhostEntries(Vec x)
{
  VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
}

// ###################################################################
/**Gets a local raw view of a ghost vector.*/
chi_math::PETScUtils::GhostVecLocalRaw
chi_math::PETScUtils::GetGhostVectorLocalViewRead(Vec x)
{
  Vec x_localized;
  VecGhostGetLocalForm(x, &x_localized);
  const double* x_localized_raw;

  VecGetArrayRead(x_localized, &x_localized_raw);

  GhostVecLocalRaw local_data;
  local_data.x_localized = x_localized;
  local_data.x_localized_raw = (double*)x_localized_raw;

  return local_data;
}

// ###################################################################
/**Gets a local raw view of a ghost vector.*/
void chi_math::PETScUtils::RestoreGhostVectorLocalViewRead(
  Vec x, GhostVecLocalRaw& local_data)
{
  VecRestoreArrayRead(local_data.x_localized,
                      (const double**)&local_data.x_localized_raw);
  VecGhostRestoreLocalForm(x, &local_data.x_localized);
}
