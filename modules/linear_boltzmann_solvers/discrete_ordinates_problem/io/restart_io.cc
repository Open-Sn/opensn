// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"

namespace opensn
{

bool
DiscreteOrdinatesProblemIO::ReadRestartData(DiscreteOrdinatesProblem& do_problem, hid_t file_id)
{
  bool success = true;

  bool file_time_dependent = do_problem.IsTimeDependent();
  if (H5Aexists(file_id, "time_dependent") > 0)
    success &= H5ReadAttribute<bool>(file_id, "time_dependent", file_time_dependent);

  OpenSnInvalidArgumentIf(file_time_dependent != do_problem.IsTimeDependent(),
                          do_problem.GetName() +
                            ": restart time-dependent mode does not match the configured "
                            "problem mode.");

  int gs_id = 0;
  for (const auto& gs : do_problem.GetGroupsets())
  {
    if (gs.angle_agg)
    {
      const auto delayed_name = "delayed_psi_old_gs" + std::to_string(gs_id);
      if (H5Has(file_id, delayed_name))
      {
        std::vector<double> psi;
        success &= H5ReadDataset1D<double>(file_id, delayed_name, psi);
        if (success)
          gs.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(psi);
      }
    }
    ++gs_id;
  }

  if (do_problem.IsTimeDependent())
  {
    const bool has_psi_old = do_problem.GetGroupsets().empty() or H5Has(file_id, "psi_old_gs0");
    OpenSnInvalidArgumentIf(not has_psi_old,
                            do_problem.GetName() +
                              ": time-dependent restart requires saved angular flux state in the "
                              "restart file.");

    auto& psi_old_local = do_problem.GetPsiOldLocal();
    auto& psi_new_local = do_problem.GetPsiNewLocal();
    psi_old_local.clear();
    psi_new_local.clear();
    psi_old_local.reserve(do_problem.GetGroupsets().size());
    psi_new_local.reserve(do_problem.GetGroupsets().size());

    for (size_t gsid = 0; gsid < do_problem.GetGroupsets().size(); ++gsid)
    {
      std::vector<double> psi_old;
      std::vector<double> psi_new;
      const auto old_name = "psi_old_gs" + std::to_string(gsid);
      const auto new_name = "psi_new_gs" + std::to_string(gsid);

      success &= H5ReadDataset1D<double>(file_id, old_name, psi_old);
      if (H5Has(file_id, new_name))
        success &= H5ReadDataset1D<double>(file_id, new_name, psi_new);
      else
        psi_new = psi_old;

      psi_old_local.push_back(std::move(psi_old));
      psi_new_local.push_back(std::move(psi_new));
    }
  }

  return success;
}

bool
DiscreteOrdinatesProblemIO::WriteRestartData(const DiscreteOrdinatesProblem& do_problem,
                                             hid_t file_id)
{
  bool success = H5CreateAttribute<bool>(file_id, "time_dependent", do_problem.IsTimeDependent());

  if (do_problem.GetOptions().write_delayed_psi_to_restart)
  {
    int gs_id = 0;
    for (const auto& gs : do_problem.GetGroupsets())
    {
      if (gs.angle_agg)
      {
        const auto psi = gs.angle_agg->GetOldDelayedAngularDOFsAsSTLVector();
        if (not psi.empty())
        {
          const auto delayed_name = "delayed_psi_old_gs" + std::to_string(gs_id);
          success &= H5WriteDataset1D<double>(file_id, delayed_name, psi);
        }
      }
      ++gs_id;
    }
  }

  if (do_problem.GetOptions().save_angular_flux)
  {
    const auto& psi_old_local = do_problem.GetPsiOldLocal();
    const auto& psi_new_local = do_problem.GetPsiNewLocal();

    for (size_t gsid = 0; gsid < psi_old_local.size(); ++gsid)
      success &= H5WriteDataset1D<double>(
        file_id, "psi_old_gs" + std::to_string(gsid), psi_old_local.at(gsid));

    for (size_t gsid = 0; gsid < psi_new_local.size(); ++gsid)
      success &= H5WriteDataset1D<double>(
        file_id, "psi_new_gs" + std::to_string(gsid), psi_new_local.at(gsid));
  }

  return success;
}

} // namespace opensn
