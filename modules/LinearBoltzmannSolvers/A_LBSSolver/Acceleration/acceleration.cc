#include "acceleration.h"

#include "mesh/SweepUtilities/SweepBoundary/sweep_boundary.h"
#include "physics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"

#include "chi_runtime.h"
#include "chi_log_exceptions.h"

namespace lbs::acceleration
{
std::map<uint64_t, BoundaryCondition>
TranslateBCs(const std::map<uint64_t, SwpBndryPtr>& sweep_boundaries,
             bool vaccum_bcs_are_dirichlet /*=true*/)
{
  typedef chi_mesh::sweep_management::BoundaryType SwpBndryType;
  typedef lbs::acceleration::BoundaryCondition BC;
  typedef lbs::acceleration::BCType BCType;

  std::map<uint64_t, BC> bcs;
  for (auto& [bid, lbs_bndry] : sweep_boundaries)
  {
    if (lbs_bndry->Type() == SwpBndryType::REFLECTING) bcs[bid] = {BCType::ROBIN, {0.0, 1.0, 0.0}};
    else if (lbs_bndry->Type() == SwpBndryType::INCIDENT_VACCUUM)
      if (vaccum_bcs_are_dirichlet) bcs[bid] = {BCType::DIRICHLET, {0.0, 0.0, 0.0}};
      else
        bcs[bid] = {BCType::ROBIN, {0.25, 0.5}};
    else // dirichlet
      bcs[bid] = {BCType::DIRICHLET, {0.0, 0.0, 0.0}};
  }

  return bcs;
}

std::map<int, Multigroup_D_and_sigR>
PackGroupsetXS(const std::map<int, MGXSPtr>& matid_to_xs_map,
               int first_grp_index,
               int last_group_index)
{
  const int num_gs_groups = last_group_index - first_grp_index + 1;
  ChiInvalidArgumentIf(num_gs_groups < 0, "last_grp_index must be >= first_grp_index");

  typedef lbs::acceleration::Multigroup_D_and_sigR MGXS;
  typedef std::map<int, lbs::acceleration::Multigroup_D_and_sigR> MatID2XSMap;
  MatID2XSMap matid_2_mgxs_map;
  for (const auto& matid_xs_pair : matid_to_xs_map)
  {
    const auto& mat_id = matid_xs_pair.first;
    const auto& xs = matid_xs_pair.second;

    std::vector<double> D(num_gs_groups, 0.0);
    std::vector<double> sigma_r(num_gs_groups, 0.0);

    size_t g = 0;
    const auto& diffusion_coeff = xs->DiffusionCoefficient();
    const auto& sigma_removal = xs->SigmaRemoval();
    for (size_t gprime = first_grp_index; gprime <= last_group_index; ++gprime)
    {
      D[g] = diffusion_coeff[gprime];
      sigma_r[g] = sigma_removal[gprime];
      ++g;
    } // for g

    matid_2_mgxs_map.insert(std::make_pair(mat_id, MGXS{D, sigma_r}));
  }

  return matid_2_mgxs_map;
}

} // namespace lbs::acceleration
