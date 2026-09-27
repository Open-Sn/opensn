// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"

namespace opensn
{
namespace
{

bool
ReadSizedDoubleVector(hid_t file_id,
                      const std::string& dataset_name,
                      std::vector<double>& destination,
                      size_t expected_size,
                      const std::string& problem_name)
{
  std::vector<double> values;
  const bool success = H5ReadDataset1D<double>(file_id, dataset_name, values);
  if (success)
  {
    OpenSnInvalidArgumentIf(values.size() != expected_size,
                            problem_name + ": restart dataset `" + dataset_name + "` has size " +
                              std::to_string(values.size()) + " but expected " +
                              std::to_string(expected_size) + ".");
    destination = std::move(values);
  }
  return success;
}

} // namespace

bool
DiscreteOrdinatesProblemIO::ReadRestartData(DiscreteOrdinatesProblem& do_problem,
                                            hid_t file_id,
                                            bool allow_transient_initialization_from_steady)
{
  bool success = true;

  bool file_time_dependent = do_problem.IsTimeDependent();
  if (H5Aexists(file_id, "time_dependent") > 0)
    success &= H5ReadAttribute<bool>(file_id, "time_dependent", file_time_dependent);

  const bool steady_to_transient_initial_condition = allow_transient_initialization_from_steady and
                                                     do_problem.IsTimeDependent() and
                                                     (not file_time_dependent);

  OpenSnInvalidArgumentIf(file_time_dependent != do_problem.IsTimeDependent() and
                            not steady_to_transient_initial_condition,
                          do_problem.GetName() +
                            ": restart time-dependent mode does not match the configured "
                            "problem mode.");

  int gs_id = 0;
  const bool full_time_dependent_restart =
    do_problem.IsTimeDependent() and not steady_to_transient_initial_condition;
  for (const auto& gs : do_problem.GetGroupsets())
  {
    if (gs.angle_agg)
    {
      const auto delayed_name = "delayed_psi_old_gs" + std::to_string(gs_id);
      const size_t expected_delayed_size = gs.angle_agg->GetNumDelayedAngularDOFs().first;
      OpenSnInvalidArgumentIf(full_time_dependent_restart and expected_delayed_size > 0 and
                                not H5Has(file_id, delayed_name),
                              do_problem.GetName() + ": time-dependent restart requires `" +
                                delayed_name + ".");
      if (H5Has(file_id, delayed_name))
      {
        std::vector<double> psi;
        success &= ReadSizedDoubleVector(
          file_id, delayed_name, psi, expected_delayed_size, do_problem.GetName());
        if (success)
          gs.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(psi);
      }
    }
    ++gs_id;
  }

  if (do_problem.IsTimeDependent())
  {
    if (not steady_to_transient_initial_condition)
    {
      auto& psi_old_local = do_problem.GetPsiOldLocal();
      auto& psi_new_local = do_problem.GetPsiNewLocal();
      std::vector<size_t> expected_psi_old_sizes;
      std::vector<size_t> expected_psi_new_sizes;
      expected_psi_old_sizes.reserve(psi_old_local.size());
      expected_psi_new_sizes.reserve(psi_new_local.size());
      for (const auto& psi : psi_old_local)
        expected_psi_old_sizes.push_back(psi.size());
      for (const auto& psi : psi_new_local)
        expected_psi_new_sizes.push_back(psi.size());
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
        const bool has_psi_old = H5Has(file_id, old_name);
        const bool has_psi_new = H5Has(file_id, new_name);
        const size_t expected_old_size = expected_psi_old_sizes.at(gsid);
        const size_t expected_new_size = expected_psi_new_sizes.at(gsid);

        OpenSnInvalidArgumentIf(
          not has_psi_old and not has_psi_new,
          do_problem.GetName() +
            ": time-dependent restart requires saved angular flux state for groupset " +
            std::to_string(gsid) + ".");

        if (has_psi_old)
          success &= ReadSizedDoubleVector(
            file_id, old_name, psi_old, expected_old_size, do_problem.GetName());
        if (has_psi_new)
          success &= ReadSizedDoubleVector(
            file_id, new_name, psi_new, expected_new_size, do_problem.GetName());
        else
          psi_new = psi_old;
        if (psi_old.empty())
          psi_old = psi_new;

        psi_old_local.push_back(std::move(psi_old));
        psi_new_local.push_back(std::move(psi_new));
      }
    }
  }

  // A steady-state restart that includes the angular flux supplies it for a transient initial
  // condition. The switch to time-dependent mode keeps a stored angular flux instead of
  // reconstructing it from the flux moments.
  if (not full_time_dependent_restart and do_problem.GetOptions().save_angular_flux)
  {
    auto& psi_old_local = do_problem.GetPsiOldLocal();
    auto& psi_new_local = do_problem.GetPsiNewLocal();
    for (size_t gsid = 0; gsid < psi_new_local.size(); ++gsid)
    {
      const auto new_name = "psi_new_gs" + std::to_string(gsid);
      if (psi_new_local[gsid].empty() or not H5Has(file_id, new_name))
        continue;
      std::vector<double> psi_new;
      success &= ReadSizedDoubleVector(
        file_id, new_name, psi_new, psi_new_local[gsid].size(), do_problem.GetName());
      if (success)
      {
        psi_new_local[gsid] = psi_new;
        if (gsid < psi_old_local.size() and psi_old_local[gsid].size() == psi_new.size())
          psi_old_local[gsid] = std::move(psi_new);
      }
    }
  }

  return success;
}

bool
DiscreteOrdinatesProblemIO::WriteRestartData(const DiscreteOrdinatesProblem& do_problem,
                                             hid_t file_id)
{
  bool success = H5CreateAttribute<bool>(file_id, "time_dependent", do_problem.IsTimeDependent());

  if (do_problem.GetOptions().restart.write_delayed_psi)
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

  if (do_problem.GetOptions().save_angular_flux and
      do_problem.GetOptions().restart.write_angular_flux)
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
