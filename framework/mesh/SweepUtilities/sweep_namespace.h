#pragma once

#include "framework/mesh/chi_mesh.h"
#include <set>

#include <memory>

namespace chi
{
class DirectedGraph;
}

// ###################################################################
namespace chi_mesh
{
namespace sweep_management
{

enum class FaceOrientation : short
{
  PARALLEL = -1,
  INCOMING = 0,
  OUTGOING = 1
};

struct STDG;     ///< Global Sweep Plane Ordering
struct SPLS;     ///< Sweep Plane Local Subgrid
class AAH_FLUDS; ///< Auxiliary Flux Data Structure
class SPDS;      ///< Sweep Plane Data Structure

struct Task
{
  unsigned int num_dependencies_;
  std::vector<uint64_t> successors_;
  uint64_t reference_id_;
  const chi_mesh::Cell* cell_ptr_;
  bool completed_ = false;
};

class AAH_ASynchronousCommunicator;
class AngleSet;
class AngleSetGroup;
class AngleAggregation;

class SweepChunk;

class SweepScheduler;

void CommunicateLocationDependencies(const std::vector<int>& location_dependencies,
                                     std::vector<std::vector<int>>& global_dependencies);

void PrintSweepOrdering(SPDS* sweep_order, MeshContinuumPtr vol_continuum);

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
typedef AngleSetStatus ExecutionPermission;
} // namespace sweep_management
} // namespace chi_mesh
