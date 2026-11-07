// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <mpi.h>
#include <string>
#include <vector>

namespace opensn
{

/**
 * Selects the GPU that this MPI rank uses and makes PETSc use the same GPU.
 *
 * Ranks on a node that see the same set of GPUs, in any order, share them round-robin in node-local
 * rank order.
 * A scheduler or launcher that gives each rank its own GPUs, by restricting the GPUs a rank can see
 * (for example with `CUDA_VISIBLE_DEVICES` or `ROCR_VISIBLE_DEVICES`), puts each rank in a set of
 * its own, so the rank uses the first GPU it was given.
 *
 * PETSc selects its device while it initializes, so when PETSc has a GPU backend this sets its
 * `-device_select_cuda` or `-device_select_hip` option, whichever matches that backend, to the
 * selected device. Selecting a different PETSc device, for example through `PETSC_OPTIONS`, makes
 * PETSc and OpenSn use different GPUs.
 *
 * Call once per process, after MPI is initialized and before `PetscInitialize` or any GPU work.
 * Collective over `comm`. Does nothing without CUDA or HIP support or on a rank that sees no GPU.
 * With SYCL, the backend selects the device.
 */
void BindDevice(MPI_Comm comm);

/**
 * Makes the GPU that BindDevice selected the current device of the calling thread. The current
 * device is per thread, so every thread other than the one that called BindDevice calls this
 * before doing GPU work. Does nothing if BindDevice selected no GPU.
 */
void UseBoundDevice();

/**
 * Returns the index, in the device order of node-local rank `node_rank`, of the GPU that the rank
 * uses, or -1 when it sees none. `visible_gpus[i]` lists the identifiers, as GetVisibleGPUIDs()
 * returns them, of the GPUs that node-local rank i sees, in its device order. Ranks that see the
 * same set of GPUs, in any order, share them round-robin in node-local rank order, taking them in
 * the order the first of those ranks lists them.
 */
int SelectDeviceIndex(const std::vector<std::vector<std::string>>& visible_gpus, int node_rank);

/**
 * Returns an identifier of each GPU this process can see, in device order. The identifier combines
 * the GPU's PCI bus ID and UUID, so it differs between the GPUs on a node and between the MIG
 * instances of one GPU, and is the same in every process that sees the GPU. Returns an empty list
 * without CUDA or HIP support, and with SYCL, whose backend selects its device itself.
 */
std::vector<std::string> GetVisibleGPUIDs();

/// Makes `device`, an index into GetVisibleGPUIDs(), the current GPU of the calling thread.
void SetCurrentGPU(int device);

} // namespace opensn
