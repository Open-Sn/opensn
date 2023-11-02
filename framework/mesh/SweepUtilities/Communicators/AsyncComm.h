#ifndef CHITECH_ASYNCCOMM_H
#define CHITECH_ASYNCCOMM_H

#include <cstddef>
#include <cstdint>
#include <vector>

namespace chi
{
class ChiMPICommunicatorSet;
}

namespace chi_mesh::sweep_management
{

class FLUDS;

class AsynchronousCommunicator
{
public:
  explicit AsynchronousCommunicator(FLUDS& fluds, const chi::ChiMPICommunicatorSet& comm_set);
  virtual ~AsynchronousCommunicator() = default;

  /**Obtains a data vector holding a spot into which outgoing data can be
   * written.*/
  virtual std::vector<double>& InitGetDownwindMessageData(int location_id,
                                                          uint64_t cell_global_id,
                                                          unsigned int face_id,
                                                          size_t angle_set_id,
                                                          size_t data_size);

protected:
  FLUDS& fluds_;
  const chi::ChiMPICommunicatorSet& comm_set_;
};

} // namespace chi_mesh::sweep_management

#endif // CHITECH_ASYNCCOMM_H
