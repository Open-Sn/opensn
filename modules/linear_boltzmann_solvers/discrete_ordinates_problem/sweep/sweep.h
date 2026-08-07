// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <set>
#include <memory>
#include <cstdint>
#include <string>

namespace opensn
{

class Cell;
class SPDS;
class MeshContinuum;

inline constexpr double FACE_ORIENTATION_TOLERANCE = 1.0e-12;

enum class FaceOrientation : short
{
  PARALLEL = -1,
  INCOMING = 0,
  OUTGOING = 1
};

enum class SweepKind
{
  AAH,
  CBC
};

/// Parses the `sweep_type` input parameter into a SweepKind.
SweepKind ParseSweepKind(const std::string& sweep_type, const std::string& problem_name);

enum class AngleSetStatus
{
  NOT_FINISHED = 0,
  FINISHED = 1,
  RECEIVING = 2,
  READY_TO_EXECUTE = 3,
  EXECUTE = 4,
  NO_EXEC_IF_READY = 5,
  MESSAGES_SENT = 6,
  MESSAGES_PENDING = 7
};

struct Task
{
  unsigned int num_dependencies;
  std::vector<std::uint32_t> successors;
  uint64_t reference_id;
  const Cell* cell_ptr;
  bool completed = false;
};

/// Stage Task Dependency Graphs
struct STDG
{
  std::vector<int> item_id;
};

/// Communicates location by location dependencies.
void CommunicateLocationDependencies(const std::vector<int>& location_dependencies,
                                     std::vector<std::vector<int>>& global_dependencies);

/// Batched variant of CommunicateLocationDependencies. Gathers dependencies for all N SPDS in 2
/// MPI collectives instead of 2*N. Returns global_deps[s][rank] for each SPDS s and MPI rank.
std::vector<std::vector<std::vector<int>>>
BatchCommunicateLocationDependencies(const std::vector<std::vector<int>>& local_deps);

/// Print a sweep ordering to file.
void PrintSweepOrdering(SPDS* sweep_order, std::shared_ptr<MeshContinuum> vol_continuum);

} // namespace opensn
