// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "caribou/caribou.h"

namespace crb = caribou;

namespace opensn
{

void
AAH_SPDS::CopySPLSDataOnDevice()
{
  // compute offset and total size to allocate on GPU
  contiguous_offset_.reserve(levelized_spls_.size());
  std::size_t total_size = 0;
  for (const std::vector<uint32_t>& level : levelized_spls_)
  {
    contiguous_offset_.push_back(total_size);
    total_size += level.size();
  }
  // allocate and copy data
  crb::DeviceMemory<std::uint32_t> device_levelized_spls(total_size);
  for (std::size_t l = 0; l < levelized_spls_.size(); ++l)
  {
    const std::vector<std::uint32_t>& level = levelized_spls_[l];
    crb::HostVector<std::uint32_t> host_level(level.begin(), level.end());
    crb::copy(device_levelized_spls, host_level, level.size(), 0, contiguous_offset_[l]);
  }
  // release ownership back to the class
  device_levelized_spls_ = device_levelized_spls.release();
}

void
AAH_SPDS::FreeDeviceData()
{
  if (device_levelized_spls_)
  {
    crb::DeviceMemory<std::uint32_t> device_levelized_spls(device_levelized_spls_);
    device_levelized_spls.reset();
    contiguous_offset_.clear();
  }
}

} // namespace opensn
