#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_curvilinear_solver/sweep_chunks/lbs_curvilinear_sweep_chunk_pwl.h"

#include "framework/math/quadratures/curvilinear_angular_quadrature.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/logging/log_exceptions.h"

namespace opensn
{
namespace lbs
{

SweepChunkPWLRZ::SweepChunkPWLRZ(
  const MeshContinuum& grid,
  const SpatialDiscretization& discretization_primary,
  const std::vector<lbs::UnitCellMatrices>& unit_cell_matrices,
  const std::vector<lbs::UnitCellMatrices>& secondary_unit_cell_matrices,
  std::vector<lbs::CellLBSView>& cell_transport_views,
  std::vector<double>& destination_phi,
  std::vector<double>& destination_psi,
  const std::vector<double>& source_moments,
  lbs::LBSGroupset& groupset,
  const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
  int num_moments,
  int max_num_cell_dofs)
  : AAH_SweepChunk(grid,
                   discretization_primary,
                   unit_cell_matrices,
                   cell_transport_views,
                   destination_phi,
                   destination_psi,
                   source_moments,
                   groupset,
                   xs,
                   num_moments,
                   max_num_cell_dofs),
    secondary_unit_cell_matrices_(secondary_unit_cell_matrices),
    unknown_manager_(),
    psi_sweep_(),
    normal_vector_boundary_()
{
  const auto curvilinear_product_quadrature =
    std::dynamic_pointer_cast<CurvilinearAngularQuadrature>(groupset_.quadrature_);

  if (curvilinear_product_quadrature == nullptr)
    throw std::invalid_argument("D_DO_RZ_SteadyState::SweepChunkPWL::SweepChunkPWL : "
                                "invalid angular quadrature");

  //  configure unknown manager for quantities that depend on polar level
  const size_t dir_map_size = curvilinear_product_quadrature->GetDirectionMap().size();
  for (size_t m = 0; m < dir_map_size; ++m)
    unknown_manager_.AddUnknown(UnknownType::VECTOR_N, groupset_.groups_.size());

  //  allocate storage for sweeping dependency
  const unsigned int n_dof = discretization_primary.GetNumLocalDOFs(unknown_manager_);
  psi_sweep_.resize(n_dof);

  //  initialise mappings from direction linear index
  for (const auto& dir_set : curvilinear_product_quadrature->GetDirectionMap())
    for (const auto& dir_idx : dir_set.second)
      map_polar_level_.emplace(dir_idx, dir_set.first);

  //  set normal vector for symmetric boundary condition
  const int d = (grid_.Attributes() & DIMENSION_1) ? 2 : 0;
  normal_vector_boundary_ = Vector3(0.0, 0.0, 0.0);
  normal_vector_boundary_(d) = 1;
}

} // namespace lbs
} // namespace opensn
