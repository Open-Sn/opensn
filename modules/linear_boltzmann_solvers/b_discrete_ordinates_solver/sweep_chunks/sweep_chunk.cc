#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/logging/log_exceptions.h"

namespace opensn
{
namespace lbs
{

SweepChunk::SweepChunk(std::vector<double>& destination_phi,
                       std::vector<double>& destination_psi,
                       const MeshContinuum& grid,
                       const SpatialDiscretization& discretization,
                       const std::vector<UnitCellMatrices>& unit_cell_matrices,
                       std::vector<lbs::CellLBSView>& cell_transport_views,
                       const std::vector<double>& source_moments,
                       const LBSGroupset& groupset,
                       const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                       int num_moments,
                       int max_num_cell_dofs)
  : opensn::SweepChunk(destination_phi, destination_psi),
    grid_(grid),
    grid_fe_view_(discretization),
    unit_cell_matrices_(unit_cell_matrices),
    grid_transport_view_(cell_transport_views),
    q_moments_(source_moments),
    groupset_(groupset),
    xs_(xs),
    num_moments_(num_moments),
    save_angular_flux_(not destination_psi.empty()),
    groupset_angle_group_stride_(groupset_.psi_uk_man_.NumberOfUnknowns() *
                                 groupset_.groups_.size()),
    groupset_group_stride_(groupset_.groups_.size()),
    max_num_cell_dofs_(max_num_cell_dofs)
{}

} // namespace lbs
} // namespace opensn
