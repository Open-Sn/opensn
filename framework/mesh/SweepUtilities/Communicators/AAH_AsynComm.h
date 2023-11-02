#ifndef CHI_AAH_ASYNCOMM_H
#define CHI_AAH_ASYNCOMM_H

#include "mesh/SweepUtilities/sweep_namespace.h"
#include "chi_mpi.h"

#include "AsyncComm.h"

typedef unsigned long long int u_ll_int;

namespace chi
{
class ChiMPICommunicatorSet;
}

namespace chi_mesh::sweep_management
{

class FLUDS;

// ###################################################################
/**Handles the swift communication of interprocess communication
 * related to sweeping.*/
class AAH_ASynchronousCommunicator : public AsynchronousCommunicator
{
private:
  const size_t num_groups_;
  const size_t num_angles_;

  bool done_sending;
  bool data_initialized;
  bool upstream_data_initialized;

  u_ll_int EAGER_LIMIT = 32000;

  std::vector<int> prelocI_message_count;
  std::vector<int> deplocI_message_count;
  std::vector<int> delayed_prelocI_message_count;

  std::vector<std::vector<u_ll_int>> prelocI_message_size;
  std::vector<std::vector<u_ll_int>> deplocI_message_size;
  std::vector<std::vector<u_ll_int>> delayed_prelocI_message_size;

  std::vector<std::vector<u_ll_int>> prelocI_message_blockpos;
  std::vector<std::vector<u_ll_int>> deplocI_message_blockpos;
  std::vector<std::vector<u_ll_int>> delayed_prelocI_message_blockpos;

  std::vector<std::vector<bool>> prelocI_message_received;
  std::vector<std::vector<bool>> delayed_prelocI_message_received;

  std::vector<std::vector<MPI_Request>> deplocI_message_request;

public:
  int max_num_mess;

  AAH_ASynchronousCommunicator(FLUDS& fluds,
                               size_t num_groups,
                               size_t num_angles,
                               int sweep_eager_limit,
                               const chi::ChiMPICommunicatorSet& in_comm_set);
  bool DoneSending() const;
  void InitializeDelayedUpstreamData();
  void InitializeLocalAndDownstreamBuffers();
  void SendDownstreamPsi(int angle_set_num);
  bool ReceiveDelayedData(int angle_set_num);
  void ClearDownstreamBuffers();
  AngleSetStatus ReceiveUpstreamPsi(int angle_set_num);
  void ClearLocalAndReceiveBuffers();
  void Reset();

protected:
  void BuildMessageStructure();
};
} // namespace chi_mesh::sweep_management
#endif // CHI_AAH_ASYNCOMM_H