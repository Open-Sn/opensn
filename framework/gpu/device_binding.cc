// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/gpu/device_binding.h"
#include "petscsys.h"
#include <algorithm>
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
SelectDeviceIndex(const std::vector<std::vector<std::string>>& visible_gpus, int node_rank)
{
  const auto& own = visible_gpus.at(node_rank);
  if (own.empty())
    return -1;

  // Ranks see the same set of GPUs when their sorted identifier lists are equal. They take the GPUs
  // in the order the first of them lists them, so ranks that list them in different orders still
  // select different GPUs.
  const auto sorted = [](std::vector<std::string> ids)
  {
    std::sort(ids.begin(), ids.end());
    return ids;
  };
  const auto own_set = sorted(own);
  const std::vector<std::string>* reference = &own;
  std::size_t position = 0;
  for (int i = 0; i < node_rank; ++i)
    if (sorted(visible_gpus[i]) == own_set)
    {
      if (position == 0)
        reference = &visible_gpus[i];
      ++position;
    }

  const auto& selected = (*reference)[position % own.size()];
  return static_cast<int>(std::find(own.begin(), own.end(), selected) - own.begin());
}

void
BindDevice([[maybe_unused]] MPI_Comm comm)
{
#ifdef __OPENSN_WITH_GPU__
  // The GPUs this rank sees
  const auto ids = GetVisibleGPUIDs();
  std::string device_set;
  for (const auto& id : ids)
  {
    device_set += id;
    device_set += ';';
  }

  // Gather the GPUs seen by every rank on this node
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

  // Each identifier is terminated by ';'
  std::vector<std::vector<std::string>> visible_gpus(node_size);
  for (int i = 0; i < node_size; ++i)
  {
    auto begin = static_cast<std::size_t>(offsets[i]);
    const auto end = begin + static_cast<std::size_t>(lengths[i]);
    while (begin < end)
    {
      const auto separator = all_sets.find(';', begin);
      visible_gpus[i].push_back(all_sets.substr(begin, separator - begin));
      begin = separator + 1;
    }
  }

  const int device = SelectDeviceIndex(visible_gpus, node_rank);
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
GetVisibleGPUIDs()
{
  return {};
}

void
SetCurrentGPU(int /*device*/)
{
}
#endif

} // namespace opensn
