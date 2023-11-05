#include "opensn/framework/physics/PhysicsMaterial/MultiGroupXS/adjoint_mgxs.h"

chi_physics::AdjointMGXS::AdjointMGXS(const MultiGroupXS& xs) : xs_(xs)
{
  // transpose transfer matrices
  for (unsigned int ell = 0; ell <= xs_.ScatteringOrder(); ++ell)
  {
    const auto& S_ell = xs_.TransferMatrix(ell);
    chi_math::SparseMatrix S_ell_transpose(xs_.NumGroups(), xs_.NumGroups());
    for (size_t g = 0; g < xs_.NumGroups(); ++g)
    {
      const size_t row_len = S_ell.rowI_indices_[g].size();
      const size_t* col_ptr = S_ell.rowI_indices_[g].data();
      const double* val_ptr = S_ell.rowI_values_[g].data();

      for (size_t j = 0; j < row_len; ++j)
        S_ell_transpose.Insert(*col_ptr++, g, *val_ptr++);
    }
    transposed_transfer_matrices_.push_back(S_ell_transpose);
  } // for ell

  // transpose production matrices
  if (xs_.IsFissionable())
  {
    const auto& F = xs_.ProductionMatrix();
    for (size_t g = 0; g < xs_.NumGroups(); ++g)
    {
      std::vector<double> F_g_transpose;
      for (size_t gp = 0; gp < xs_.NumGroups(); ++gp)
        F_g_transpose.emplace_back(F[gp][g]);
      transposed_production_matrices_.push_back(F_g_transpose);
    }
  }
}
