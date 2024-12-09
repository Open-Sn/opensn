#include "framework/math/parallel_vector/parallel_stl_vector.h"
#include "framework/math/parallel_vector/ghosted_parallel_stl_vector.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/lib/console.h"
#include "LuaBridge/LuaBridge.h"

using namespace opensn;

namespace unit_tests
{

void
math_Test02_ParallelVector()
{
  OpenSnLogicalErrorIf(opensn::mpi_comm.size() != 2, "Requires 2 processors");

  opensn::log.Log() << "Testing ParallelSTLVector" << std::endl;

  ParallelSTLVector vec(5, 10, opensn::mpi_comm);

  if (opensn::mpi_comm.rank() == 0)
    vec.SetValue(5, 2.0, VecOpType::SET_VALUE);
  else
    vec.SetValue(0, 1.0, VecOpType::SET_VALUE);
  vec.Assemble();

  opensn::log.LogAll() << "vec after assembly: " << vec.PrintStr() << std::endl;

  opensn::log.Log() << "Testing GhostedParallelSTLVector" << std::endl;

  const int64_t ghost_id = opensn::mpi_comm.rank() == 0 ? 5 : 4;
  GhostedParallelSTLVector ghost_vec(5, 10, {ghost_id}, opensn::mpi_comm);

  opensn::log.LogAll() << "Number of Ghosts: " << ghost_vec.GetNumGhosts() << std::endl;

  if (opensn::mpi_comm.rank() == 0)
    ghost_vec.SetValue(5, 2.0, VecOpType::SET_VALUE);
  else
    ghost_vec.SetValue(4, 1.0, VecOpType::SET_VALUE);
  ghost_vec.Assemble();
  ghost_vec.CommunicateGhostEntries();

  opensn::log.LogAll() << "Ghost vec after communicate: " << ghost_vec.PrintStr() << std::endl;

  {
    std::stringstream outstr;
    // const auto& raw_vals = ghost_vec.MakeLocalVector();
    const double* data = ghost_vec.Data();
    for (size_t i = 0; i < ghost_vec.GetLocalSizeWithGhosts(); ++i)
      outstr << data[i] << " ";
    opensn::log.LogAll() << "Ghost vec raw values: " << outstr.str();
  }

  {
    std::stringstream outstr;
    const auto made_vals = ghost_vec.MakeLocalVector();
    for (double val : made_vals)
      outstr << val << " ";
    opensn::log.LogAll() << "Ghost vec make-local values: " << outstr.str();
  }

  opensn::log.LogAll() << "Parallel vector norm-1: " << vec.ComputeNorm(opensn::NormType::L1_NORM);
  opensn::log.LogAll() << "Parallel vector norm-2: " << vec.ComputeNorm(opensn::NormType::L2_NORM);
  opensn::log.LogAll() << "Parallel vector norm-inf: "
                       << vec.ComputeNorm(opensn::NormType::LINF_NORM);

  opensn::log.LogAll() << "Ghost vector norm-1: "
                       << ghost_vec.ComputeNorm(opensn::NormType::L1_NORM);
  opensn::log.LogAll() << "Ghost vector norm-2: "
                       << ghost_vec.ComputeNorm(opensn::NormType::L2_NORM);
  opensn::log.LogAll() << "Ghost vector norm-inf: "
                       << ghost_vec.ComputeNorm(opensn::NormType::LINF_NORM);

  opensn::log.Log() << "Testing ParallelSTLVector "
                    << "ADD_VALUE and CopyValues" << std::endl;
  ParallelSTLVector vec2(5, 10, opensn::mpi_comm);

  vec2.CopyLocalValues(vec);

  if (opensn::mpi_comm.rank() == 0)
    vec2.SetValue(5, 2.0, VecOpType::ADD_VALUE);
  else
    vec2.SetValue(0, 1.0, VecOpType::ADD_VALUE);
  vec2.Assemble();

  opensn::log.LogAll() << "vec2 after assembly: " << vec2.PrintStr() << std::endl;

  opensn::log.Log() << "Testing ParallelSTLVector "
                    << "SetValues" << std::endl;
  ParallelSTLVector vec3(5, 10, opensn::mpi_comm);

  if (opensn::mpi_comm.rank() == 0)
    vec3.SetValues({5, 6}, {2.0, 3.0}, VecOpType::ADD_VALUE);
  else
    vec3.SetValues({0, 1}, {1.0, 4.0}, VecOpType::ADD_VALUE);
  vec3.Assemble();

  opensn::log.LogAll() << "vec3 after assembly: " << vec3.PrintStr() << std::endl;

  opensn::log.Log() << "Testing GhostedParallelSTLVector "
                    << "Constructed from VectorGhostCommunicator and "
                       "other utilities"
                    << std::endl;

  std::vector<int64_t> ghost_ids;
  if (opensn::mpi_comm.rank() == 0)
    ghost_ids = {5, 6};
  else
    ghost_ids = {0, 1, 3};
  VectorGhostCommunicator vgc(5, 10, ghost_ids, opensn::mpi_comm);

  GhostedParallelSTLVector ghost_vec2(vgc);

  opensn::log.LogAll() << "ghost_vec2 local size with ghosts "
                       << ghost_vec2.GetLocalSizeWithGhosts() << std::endl;

  if (opensn::mpi_comm.rank() == 0)
    ghost_vec2.SetValues({5, 6}, {6.0, 7.0}, VecOpType::ADD_VALUE);
  else
    ghost_vec2.SetValues({0, 1, 3}, {1.0, 2.0, 4.0}, VecOpType::ADD_VALUE);

  ghost_vec2.Assemble();
  ghost_vec2.CommunicateGhostEntries();

  opensn::log.LogAll() << "ghost_vec2 after assembly: " << ghost_vec2.PrintStr() << std::endl;

  {
    std::stringstream outstr;
    for (int64_t gid : ghost_vec2.GetGhostIndices())
      outstr << gid << " ";
    opensn::log.LogAll() << "ghost_vec2 ghost ids: " << outstr.str() << std::endl;
  }

  if (opensn::mpi_comm.rank() == 0)
    opensn::log.LogAll() << "ghost_vec2 mapghostA: " << ghost_vec2.MapGhostToLocal(6) << std::endl;
  else
    opensn::log.LogAll() << "ghost_vec2 mapghostA: " << ghost_vec2.MapGhostToLocal(1) << std::endl;

  {
    const auto ghosted_local = ghost_vec2.MakeGhostedLocalVector();
    std::stringstream outstr;
    for (double gid : ghosted_local)
      outstr << gid << " ";
    opensn::log.LogAll() << "ghost_vec2 MakeGhostedLocalVector: " << outstr.str() << std::endl;
  }

  if (opensn::mpi_comm.rank() == 0)
    opensn::log.LogAll() << "ghost_vec2 GetGlobalValue(local): " << ghost_vec2.GetGlobalValue(3)
                         << std::endl;
  else
    opensn::log.LogAll() << "ghost_vec2 GetGlobalValue(local): " << ghost_vec2.GetGlobalValue(6)
                         << std::endl;

  if (opensn::mpi_comm.rank() == 0)
    opensn::log.LogAll() << "ghost_vec2 GetGlobalValue(ghost): " << ghost_vec2.GetGlobalValue(6)
                         << std::endl;
  else
    opensn::log.LogAll() << "ghost_vec2 GetGlobalValue(ghost): " << ghost_vec2.GetGlobalValue(1)
                         << std::endl;
}

BIND_FUNCTION(unit_tests, math_Test02_ParallelVector);

} //  namespace unit_tests
