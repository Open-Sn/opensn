// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/lbs.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"

using namespace opensn;

namespace opensnlua
{

std::vector<std::shared_ptr<FieldFunctionGridBased>>
LBSGetScalarFieldFunctionList(std::shared_ptr<LBSProblem> lbs_problem)
{
  // Building table of handles
  std::vector<std::shared_ptr<FieldFunctionGridBased>> fflist;
  // Flux moments first
  for (int g = 0; g < lbs_problem->GetNumGroups(); g++)
  {
    for (int m = 0; m < lbs_problem->GetNumMoments(); m++)
    {
      auto ff = lbs_problem->MapPhiFieldFunction(g, m);
      auto local_ff = lbs_problem->GetFieldFunctions()[ff];

      if (m != 0)
        continue;

      fflist.push_back(local_ff);
    }
  }
  return fflist;
}

std::map<std::string, std::vector<double>>
LBSComputeLeakage(std::shared_ptr<opensn::DiscreteOrdinatesProblem> solver,
                  const std::vector<std::string>& bnd_names)
{
  // Get the supported boundaries
  const auto supported_boundary_names = opensn::DiscreteOrdinatesProblem::supported_boundary_names;
  const auto supported_boundary_ids = opensn::DiscreteOrdinatesProblem::supported_boundary_ids;

  // Get the boundaries to parse
  std::vector<uint64_t> bndry_ids;
  if (bnd_names.size() > 1)
  {
    for (auto& name : bnd_names)
      bndry_ids.push_back(supported_boundary_names.at(name));
  }
  else
    bndry_ids = solver->GetGrid()->GetUniqueBoundaryIDs();

  // Compute the leakage
  const auto leakage = solver->ComputeLeakage(bndry_ids);

  std::map<std::string, std::vector<double>> ret_val;
  for (const auto& [bid, vals] : leakage)
  {
    auto bnd_name = supported_boundary_ids.at(bid);
    ret_val.insert(std::pair<std::string, std::vector<double>>(bnd_name, vals));
  }
  return ret_val;
}

void
LBSWriteFluxMoments(std::shared_ptr<LBSProblem> lbs_problem, const std::string& file_base)
{
  LBSSolverIO::WriteFluxMoments(*lbs_problem, file_base);
}

void
LBSCreateAndWriteSourceMoments(std::shared_ptr<LBSProblem> lbs_problem,
                               const std::string& file_base)
{
  auto source_moments = lbs_problem->MakeSourceMomentsFromPhi();
  LBSSolverIO::WriteFluxMoments(*lbs_problem, file_base, source_moments);
}

void
LBSReadFluxMomentsAndMakeSourceMoments(std::shared_ptr<LBSProblem> lbs_problem,
                                       const std::string& file_base,
                                       bool single_file_flag)
{
  LBSSolverIO::ReadFluxMoments(
    *lbs_problem, file_base, single_file_flag, lbs_problem->GetExtSrcMomentsLocal());

  opensn::log.Log() << "Making source moments from flux file.";
  auto temp_phi = lbs_problem->GetPhiOldLocal();
  lbs_problem->GetPhiOldLocal() = lbs_problem->GetExtSrcMomentsLocal();
  lbs_problem->GetExtSrcMomentsLocal() = lbs_problem->MakeSourceMomentsFromPhi();
  lbs_problem->GetPhiOldLocal() = temp_phi;
}

void
LBSReadSourceMoments(std::shared_ptr<LBSProblem> lbs_problem,
                     const std::string& file_base,
                     bool single_file_flag)
{
  LBSSolverIO::ReadFluxMoments(
    *lbs_problem, file_base, single_file_flag, lbs_problem->GetExtSrcMomentsLocal());
}

void
LBSReadFluxMoments(std::shared_ptr<LBSProblem> lbs_problem,
                   const std::string& file_base,
                   bool single_file_flag)
{
  LBSSolverIO::ReadFluxMoments(*lbs_problem, file_base, single_file_flag);
}

void
LBSWriteAngularFluxes(std::shared_ptr<LBSProblem> lbs_problem, const std::string& file_base)
{
  LBSSolverIO::WriteAngularFluxes(*lbs_problem, file_base);
}

void
LBSReadAngularFluxes(std::shared_ptr<LBSProblem> lbs_problem, const std::string& file_base)
{
  LBSSolverIO::ReadAngularFluxes(*lbs_problem, file_base);
}

} // namespace opensnlua
