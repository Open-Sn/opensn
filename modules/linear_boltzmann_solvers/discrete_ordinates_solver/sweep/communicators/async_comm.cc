#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/communicators/async_comm.h"
#include "framework/logging/log_exceptions.h"

namespace opensn
{
namespace lbs
{

AsynchronousCommunicator::AsynchronousCommunicator(FLUDS& fluds, const MPICommunicatorSet& comm_set)
  : fluds_(fluds), comm_set_(comm_set)
{
}

std::vector<double>&
AsynchronousCommunicator::InitGetDownwindMessageData(int location_id,
                                                     uint64_t cell_global_id,
                                                     unsigned int face_id,
                                                     size_t angle_set_id,
                                                     size_t data_size)
{
  OpenSnLogicalError("Method not implemented");
}

} // namespace lbs
} // namespace opensn
