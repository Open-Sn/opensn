// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/physics/physics_material/multi_group_xs/adjoint_mgxs.h"

namespace opensn
{

AdjointMGXS::AdjointMGXS(const MGXS& xs) : xs_(xs)
{
  // Transpose transfer matrices
  for (unsigned int ell = 0; ell <= xs_.ScatteringOrder(); ++ell)
  {
    const auto& S_ell = xs_.TransferMatrix(ell);
    SparseMatrix S_ell_transpose(xs_.NumGroups(), xs_.NumGroups());
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

  // Transpose production matrices
  if (xs_.IsFissionable())
  {
    transposed_production_matrices_.clear();
    transposed_production_matrices_.resize(xs_.NumGroups());
    const auto& F = xs_.ProductionMatrix();
    for (size_t g = 0; g < xs_.NumGroups(); ++g)
      for (size_t gp = 0; gp < xs_.NumGroups(); ++gp)
        transposed_production_matrices_[g].push_back(F[gp][g]);
  }
}

} // namespace opensn
