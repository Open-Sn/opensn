#include "lbs_solver.h"

//###################################################################
/**Initializes common groupset items.*/
void lbs::LBSSolver::InitializeGroupsets()
{
  for (auto& groupset : groupsets_)
  {
    //================================================== Build groupset angular
    //                                                   flux unknown manager
    groupset.psi_uk_man_.unknowns_.clear();
    size_t num_angles = groupset.quadrature_->abscissae_.size();
    size_t gs_num_groups = groupset.groups_.size();
    auto& grpset_psi_uk_man = groupset.psi_uk_man_;

    const auto VarVecN = chi_math::UnknownType::VECTOR_N;
    for (unsigned int n=0; n<num_angles; ++n)
      grpset_psi_uk_man.AddUnknown(VarVecN, gs_num_groups);

    groupset.BuildDiscMomOperator(options_.scattering_order,
                                  options_.geometry_type);
    groupset.BuildMomDiscOperator(options_.scattering_order,
                                  options_.geometry_type);
    groupset.BuildSubsets();
  }//for groupset
}