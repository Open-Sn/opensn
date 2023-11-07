#include "framework/math/ParallelVector/ParallelSTLVector.h"
#include "framework/math/ParallelVector/GhostedParallelSTLVector.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "framework/console/chi_console.h"

namespace chi_unit_tests
{

chi::ParameterBlock chi_math_Test02_ParallelVector(const chi::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chi_math_Test02_ParallelVector,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/chi_math_Test02_ParallelVector);

chi::ParameterBlock
chi_math_Test02_ParallelVector(const chi::InputParameters&)
{
  using namespace chi_math;

  ChiLogicalErrorIf(Chi::mpi.process_count != 2, "Requires 2 processors");

  Chi::log.Log() << "Testing chi_math::ParallelSTLVector" << std::endl;

  ParallelSTLVector vec(5, 10, Chi::mpi.comm);

  if (Chi::mpi.location_id == 0) vec.SetValue(5, 2.0, VecOpType::SET_VALUE);
  else
    vec.SetValue(0, 1.0, VecOpType::SET_VALUE);
  vec.Assemble();

  Chi::log.LogAll() << "vec after assembly: " << vec.PrintStr() << std::endl;

  Chi::log.Log() << "Testing chi_math::GhostedParallelSTLVector" << std::endl;

  const int64_t ghost_id = Chi::mpi.location_id == 0 ? 5 : 4;
  GhostedParallelSTLVector ghost_vec(5, 10, {ghost_id}, Chi::mpi.comm);

  Chi::log.LogAll() << "Number of Ghosts: " << ghost_vec.NumGhosts() << std::endl;

  if (Chi::mpi.location_id == 0) ghost_vec.SetValue(5, 2.0, VecOpType::SET_VALUE);
  else
    ghost_vec.SetValue(4, 1.0, VecOpType::SET_VALUE);
  ghost_vec.Assemble();
  ghost_vec.CommunicateGhostEntries();

  Chi::log.LogAll() << "Ghost vec after communicate: " << ghost_vec.PrintStr() << std::endl;

  {
    std::stringstream outstr;
    // const auto& raw_vals = ghost_vec.MakeLocalVector();
    const double* data = ghost_vec.Data();
    for (size_t i = 0; i < ghost_vec.LocalSizeWithGhosts(); ++i)
      outstr << data[i] << " ";
    Chi::log.LogAll() << "Ghost vec raw values: " << outstr.str();
  }

  {
    std::stringstream outstr;
    const auto made_vals = ghost_vec.MakeLocalVector();
    for (double val : made_vals)
      outstr << val << " ";
    Chi::log.LogAll() << "Ghost vec make-local values: " << outstr.str();
  }

  Chi::log.LogAll() << "Parallel vector norm-1: " << vec.ComputeNorm(chi_math::NormType::L1_NORM);
  Chi::log.LogAll() << "Parallel vector norm-2: " << vec.ComputeNorm(chi_math::NormType::L2_NORM);
  Chi::log.LogAll() << "Parallel vector norm-inf: "
                    << vec.ComputeNorm(chi_math::NormType::LINF_NORM);

  Chi::log.LogAll() << "Ghost vector norm-1: "
                    << ghost_vec.ComputeNorm(chi_math::NormType::L1_NORM);
  Chi::log.LogAll() << "Ghost vector norm-2: "
                    << ghost_vec.ComputeNorm(chi_math::NormType::L2_NORM);
  Chi::log.LogAll() << "Ghost vector norm-inf: "
                    << ghost_vec.ComputeNorm(chi_math::NormType::LINF_NORM);

  Chi::log.Log() << "Testing chi_math::ParallelSTLVector "
                 << "ADD_VALUE and CopyValues" << std::endl;
  ParallelSTLVector vec2(5, 10, Chi::mpi.comm);

  vec2.CopyLocalValues(vec);

  if (Chi::mpi.location_id == 0) vec2.SetValue(5, 2.0, VecOpType::ADD_VALUE);
  else
    vec2.SetValue(0, 1.0, VecOpType::ADD_VALUE);
  vec2.Assemble();

  Chi::log.LogAll() << "vec2 after assembly: " << vec2.PrintStr() << std::endl;

  Chi::log.Log() << "Testing chi_math::ParallelSTLVector "
                 << "SetValues" << std::endl;
  ParallelSTLVector vec3(5, 10, Chi::mpi.comm);

  if (Chi::mpi.location_id == 0) vec3.SetValues({5, 6}, {2.0, 3.0}, VecOpType::ADD_VALUE);
  else
    vec3.SetValues({0, 1}, {1.0, 4.0}, VecOpType::ADD_VALUE);
  vec3.Assemble();

  Chi::log.LogAll() << "vec3 after assembly: " << vec3.PrintStr() << std::endl;

  Chi::log.Log() << "Testing chi_math::GhostedParallelSTLVector "
                 << "Constructed from VectorGhostCommunicator and "
                    "other utilities"
                 << std::endl;

  std::vector<int64_t> ghost_ids;
  if (Chi::mpi.location_id == 0) ghost_ids = {5, 6};
  else
    ghost_ids = {0, 1, 3};
  VectorGhostCommunicator vgc(5, 10, ghost_ids, Chi::mpi.comm);

  GhostedParallelSTLVector ghost_vec2(vgc);

  Chi::log.LogAll() << "ghost_vec2 local size with ghosts " << ghost_vec2.LocalSizeWithGhosts()
                    << std::endl;

  if (Chi::mpi.location_id == 0) ghost_vec2.SetValues({5, 6}, {6.0, 7.0}, VecOpType::ADD_VALUE);
  else
    ghost_vec2.SetValues({0, 1, 3}, {1.0, 2.0, 4.0}, VecOpType::ADD_VALUE);

  ghost_vec2.Assemble();
  ghost_vec2.CommunicateGhostEntries();

  Chi::log.LogAll() << "ghost_vec2 after assembly: " << ghost_vec2.PrintStr() << std::endl;

  {
    std::stringstream outstr;
    for (int64_t gid : ghost_vec2.GhostIndices())
      outstr << gid << " ";
    Chi::log.LogAll() << "ghost_vec2 ghost ids: " << outstr.str() << std::endl;
  }

  if (Chi::mpi.location_id == 0)
    Chi::log.LogAll() << "ghost_vec2 mapghostA: " << ghost_vec2.MapGhostToLocal(6) << std::endl;
  else
    Chi::log.LogAll() << "ghost_vec2 mapghostA: " << ghost_vec2.MapGhostToLocal(1) << std::endl;

  {
    const auto ghosted_local = ghost_vec2.MakeGhostedLocalVector();
    std::stringstream outstr;
    for (double gid : ghosted_local)
      outstr << gid << " ";
    Chi::log.LogAll() << "ghost_vec2 MakeGhostedLocalVector: " << outstr.str() << std::endl;
  }

  if (Chi::mpi.location_id == 0)
    Chi::log.LogAll() << "ghost_vec2 GetGlobalValue(local): " << ghost_vec2.GetGlobalValue(3)
                      << std::endl;
  else
    Chi::log.LogAll() << "ghost_vec2 GetGlobalValue(local): " << ghost_vec2.GetGlobalValue(6)
                      << std::endl;

  if (Chi::mpi.location_id == 0)
    Chi::log.LogAll() << "ghost_vec2 GetGlobalValue(ghost): " << ghost_vec2.GetGlobalValue(6)
                      << std::endl;
  else
    Chi::log.LogAll() << "ghost_vec2 GetGlobalValue(ghost): " << ghost_vec2.GetGlobalValue(1)
                      << std::endl;

  return chi::ParameterBlock();
}

} // namespace chi_unit_tests
