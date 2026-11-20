// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/gpu/device_binding.h"
#include "petscsys.h"
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
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

  // A GPU that a rank lists more than once, as SYCL lists a GPU that two of its backends reach,
  // counts once.
  const auto distinct = [](const std::vector<std::string>& ids)
  {
    std::vector<std::string> gpus;
    for (const auto& id : ids)
      if (std::find(gpus.begin(), gpus.end(), id) == gpus.end())
        gpus.push_back(id);
    return gpus;
  };
  // Ranks see the same set of GPUs when their sorted lists of distinct identifiers are equal. They
  // take the GPUs in the order the first of them lists them, so ranks that list them in different
  // orders still select different GPUs.
  const auto sorted = [](std::vector<std::string> ids)
  {
    std::sort(ids.begin(), ids.end());
    return ids;
  };
  auto reference = distinct(own);
  const auto own_set = sorted(reference);
  std::size_t position = 0;
  for (int i = 0; i < node_rank; ++i)
  {
    auto gpus = distinct(visible_gpus[i]);
    if (sorted(gpus) == own_set)
    {
      if (position == 0)
        reference = std::move(gpus);
      ++position;
    }
  }

  // The rank uses its first device with the selected identifier
  const auto& selected = reference[position % reference.size()];
  return static_cast<int>(std::find(own.begin(), own.end(), selected) - own.begin());
}

void
BindDevice([[maybe_unused]] MPI_Comm comm)
{
#ifdef __OPENSN_WITH_GPU__
  // A rank whose GPU query or selection fails still takes part in every collective below, so that
  // every rank reports the failure rather than waiting for the failed rank.
  std::string error;

  // The GPUs this rank sees
  std::vector<std::string> ids;
  try
  {
    ids = GetVisibleGPUIDs();
  }
  catch (const std::exception& e)
  {
    error = e.what();
  }
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

  if (error.empty())
  {
    try
    {
      const int device = SelectDeviceIndex(visible_gpus, node_rank);
      if (device >= 0)
      {
        SetCurrentGPU(device);
#if defined(PETSC_HAVE_CUDA) || defined(PETSC_HAVE_HIP)
        // By default PETSc selects (global rank) % (number of GPUs) while it initializes. PETSc
        // and OpenSn share the current device of the CUDA or HIP runtime, so they must select the
        // same one.
        if (GPUBackendIsCUDAOrHIP())
        {
#if defined(PETSC_HAVE_CUDA)
          const std::string option = "-device_select_cuda";
#else
          const std::string option = "-device_select_hip";
#endif
          if (PetscOptionsSetValue(nullptr, option.c_str(), std::to_string(device).c_str()) !=
              PETSC_SUCCESS)
            throw std::runtime_error("failed to set PETSc option " + option);
        }
#endif
        bound_device = device;
      }
    }
    catch (const std::exception& e)
    {
      error = e.what();
    }
  }

  // Every rank fails if any rank failed
  const int failed = error.empty() ? 0 : 1;
  int any_failed = 0;
  MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, comm);
  if (failed != 0)
    throw std::runtime_error("BindDevice: GPU selection failed on rank " + std::to_string(rank) +
                             ": " + error);
  if (any_failed != 0)
    throw std::runtime_error("BindDevice: GPU selection failed on another rank");
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

bool
GPUBackendIsCUDAOrHIP()
{
  return false;
}
#endif

} // namespace opensn
