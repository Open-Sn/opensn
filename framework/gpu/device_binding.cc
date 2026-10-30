// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/gpu/device_binding.h"
#include "petscsys.h"
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace opensn
{

namespace
{

// Set once by BindDevice, before any other thread starts
int bound_device = -1;

} // namespace

int
SelectDeviceIndex(const std::vector<std::string>& device_sets, int node_rank, int num_visible)
{
  if (num_visible <= 0)
    return -1;
  // Position of this rank among the node-local ranks that see the same GPUs
  const auto& own_set = device_sets.at(node_rank);
  int position = 0;
  for (int i = 0; i < node_rank; ++i)
    if (device_sets[i] == own_set)
      ++position;
  return position % num_visible;
}

void
BindDevice([[maybe_unused]] MPI_Comm comm)
{
#ifdef __OPENSN_WITH_GPU__
  // The GPUs this rank sees, identified by PCI bus ID
  const auto bus_ids = GetVisibleGPUBusIDs();
  const auto num_visible = static_cast<int>(bus_ids.size());
  std::string device_set;
  for (const auto& bus_id : bus_ids)
  {
    device_set += bus_id;
    device_set += ';';
  }

  // Gather the sets seen by every rank on this node
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm node_comm = MPI_COMM_NULL;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node_comm);
  int node_rank = 0;
  int node_size = 0;
  MPI_Comm_rank(node_comm, &node_rank);
  MPI_Comm_size(node_comm, &node_size);

  const int length = static_cast<int>(device_set.size());
  std::vector<int> lengths(node_size, 0);
  MPI_Allgather(&length, 1, MPI_INT, lengths.data(), 1, MPI_INT, node_comm);
  std::vector<int> offsets(node_size, 0);
  for (int i = 1; i < node_size; ++i)
    offsets[i] = offsets[i - 1] + lengths[i - 1];
  std::string all_sets(static_cast<std::size_t>(offsets.back() + lengths.back()), '\0');
  MPI_Allgatherv(device_set.data(),
                 length,
                 MPI_CHAR,
                 all_sets.data(),
                 lengths.data(),
                 offsets.data(),
                 MPI_CHAR,
                 node_comm);
  MPI_Comm_free(&node_comm);

  std::vector<std::string> device_sets;
  device_sets.reserve(node_size);
  for (int i = 0; i < node_size; ++i)
    device_sets.push_back(all_sets.substr(offsets[i], lengths[i]));

  const int device = SelectDeviceIndex(device_sets, node_rank, num_visible);
  if (device < 0)
    return;
  SetCurrentGPU(device);
  bound_device = device;

#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
  // By default PETSc selects (global rank) % (number of GPUs) while it initializes.
#if defined(PETSC_HAVE_CUDA)
  const std::string option = "-device_select_cuda";
#else
  const std::string option = "-device_select_hip";
#endif
  if (PetscOptionsSetValue(nullptr, option.c_str(), std::to_string(device).c_str()) !=
      PETSC_SUCCESS)
    throw std::runtime_error("BindDevice: failed to set PETSc option " + option);
#endif
#endif
}

void
UseBoundDevice()
{
  if (bound_device >= 0)
    SetCurrentGPU(bound_device);
}

#ifndef __OPENSN_WITH_GPU__
std::vector<std::string>
GetVisibleGPUBusIDs()
{
  return {};
}

void
SetCurrentGPU(int /*device*/)
{
}
#endif

} // namespace opensn
