#include "framework/mesh/sweep_utilities/communicators/async_comm.h"

#include "framework/logging/log_exceptions.h"

namespace chi_mesh::sweep_management
{

AsynchronousCommunicator::AsynchronousCommunicator(FLUDS& fluds,
                                                   const chi::ChiMPICommunicatorSet& comm_set)
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
  ChiLogicalError("Method not implemented");
}

} // namespace chi_mesh::sweep_management
