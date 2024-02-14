#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"

namespace opensn
{
namespace lbs
{

SweepChunk::SweepChunk(std::vector<double>& destination_phi,
                       std::vector<double>& destination_psi,
                       const MeshContinuum& grid,
                       const SpatialDiscretization& grid_fe_view,
                       const std::vector<UnitCellMatrices>& unit_cell_matrices,
                       std::vector<lbs::CellLBSView>& grid_transport_view,
                       const std::vector<double>& source_moments,
                       const LBSGroupset& groupset,
                       const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                       int num_moments,
                       int max_num_cell_dofs)
  : opensn::SweepChunk(destination_phi, destination_psi),
    grid_(grid),
    grid_fe_view_(grid_fe_view),
    unit_cell_matrices_(unit_cell_matrices),
    grid_transport_view_(grid_transport_view),
    source_moments_(source_moments),
    groupset_(groupset),
    xs_(xs),
    num_moments_(num_moments),
    max_num_cell_dofs_(max_num_cell_dofs),
    save_angular_flux_(not destination_psi.empty()),
    groupset_angle_group_stride_(groupset_.psi_uk_man_.NumberOfUnknowns()*groupset_.groups_.size()),
    groupset_group_stride_(groupset_.groups_.size())
{}

} // namespace lbs
} // namespace opensn
