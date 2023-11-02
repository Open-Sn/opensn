#include "petsc_utils.h"

#include "chi_log.h"

//###################################################################
/**Creates a general square matrix.
 *
 * This is a macro for:
\code
Mat A;
MatCreate(PETSC_COMM_WORLD,&A);
MatSetType(A,MATMPIAIJ);
MatSetSizes(A,local_size, local_size,
              global_size, global_size);

MatMPIAIJSetPreallocation(A,1, nullptr,
                            0, nullptr);
MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

return A;
\endcode

*/
Mat chi_math::PETScUtils::CreateSquareMatrix(int64_t local_size, int64_t global_size)
{
  Mat A;
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetType(A,MATMPIAIJ);
  MatSetSizes(A,local_size, local_size,
              global_size, global_size);

  MatMPIAIJSetPreallocation(A,1, nullptr,
                            0, nullptr);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

  return A;
}

//###################################################################
/**Creates a general square matrix.
 *
 * This is a macro for:
\code
MatCreate(PETSC_COMM_WORLD,&A);
MatSetType(A,MATMPIAIJ);
MatSetSizes(A,local_size, local_size,
              global_size, global_size);

MatMPIAIJSetPreallocation(A,1, nullptr,
                          0, nullptr);
MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
\endcode

*/
void chi_math::PETScUtils::
CreateSquareMatrix(Mat& A, int64_t local_size, int64_t global_size)
{
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetType(A,MATMPIAIJ);
  MatSetSizes(A,local_size, local_size,
              global_size, global_size);

  MatMPIAIJSetPreallocation(A,1, nullptr,
                            0, nullptr);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
}

//###################################################################
/**Initializes the sparsity pattern of a matrix.

This is a macro for:
\code
MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
                            0,nodal_nnz_off_diag.data());
MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
MatSetUp(A);
\endcode
*/
void chi_math::PETScUtils::InitMatrixSparsity(
  Mat &A,
  const std::vector<int64_t>& nodal_nnz_in_diag,
  const std::vector<int64_t>& nodal_nnz_off_diag)
{
  MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
                              0,nodal_nnz_off_diag.data());
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  MatSetUp(A);
}

//###################################################################
/**Initializes the sparsity pattern of a matrix.

This is a macro for:
\code
MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
                            0,nodal_nnz_off_diag.data());
MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
MatSetUp(A);
\endcode
*/
void chi_math::PETScUtils::InitMatrixSparsity(
  Mat &A,
  int64_t nodal_nnz_in_diag,
  int64_t nodal_nnz_off_diag)
{
  MatMPIAIJSetPreallocation(A,nodal_nnz_in_diag, nullptr,
                              nodal_nnz_off_diag, nullptr);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
}