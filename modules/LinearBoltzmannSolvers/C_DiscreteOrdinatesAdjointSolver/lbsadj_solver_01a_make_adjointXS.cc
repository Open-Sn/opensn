#include "lbsadj_solver.h"

#include "physics/PhysicsMaterial/MultiGroupXS/adjoint_mgxs.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

namespace lbs
{

void
DiscreteOrdinatesAdjointSolver::MakeAdjointXSs()
{
  //============================================= Create adjoint cross sections
  using AdjXS = chi_physics::AdjointMGXS;

  // define the actual cross-sections
  std::map<int, XSPtr> matid_to_adj_xs_map;
  for (const auto& matid_xs_pair : matid_to_xs_map_)
  {
    const auto matid = matid_xs_pair.first;
    const auto fwd_xs = std::dynamic_pointer_cast<chi_physics::MultiGroupXS>(matid_xs_pair.second);
    matid_to_adj_xs_map[matid] = std::make_shared<AdjXS>(*fwd_xs);
  } // for each mat
  matid_to_xs_map_ = std::move(matid_to_adj_xs_map);

  // reassign transport view to adjoint cross-sections
  if (grid_ptr_->local_cells.size() == cell_transport_views_.size())
    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& xs_ptr = matid_to_xs_map_[cell.material_id_];
      auto& transport_view = cell_transport_views_[cell.local_id_];

      transport_view.ReassingXS(*xs_ptr);
    }
}

} // namespace lbs