// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"

namespace opensnlua
{

std::vector<std::shared_ptr<opensn::FieldFunctionGridBased>>
LBSGetScalarFieldFunctionList(std::shared_ptr<opensn::LBSProblem> lbs_problem);

std::map<std::string, std::vector<double>>
LBSComputeLeakage(std::shared_ptr<opensn::DiscreteOrdinatesProblem> solver,
                  const std::vector<std::string>& bnd_names);

void LBSWriteFluxMoments(std::shared_ptr<opensn::LBSProblem> lbs_problem,
                         const std::string& file_base);

void LBSCreateAndWriteSourceMoments(std::shared_ptr<opensn::LBSProblem> lbs_problem,
                                    const std::string& file_base);

void LBSReadFluxMomentsAndMakeSourceMoments(std::shared_ptr<opensn::LBSProblem> lbs_problem,
                                            const std::string& file_base,
                                            bool single_file_flag);

void LBSReadSourceMoments(std::shared_ptr<opensn::LBSProblem> lbs_problem,
                          const std::string& file_base,
                          bool single_file_flag);

void LBSReadFluxMoments(std::shared_ptr<opensn::LBSProblem> lbs_problem,
                        const std::string& file_base,
                        bool single_file_flag);

void LBSWriteAngularFluxes(std::shared_ptr<opensn::LBSProblem> lbs_problem,
                           const std::string& file_base);

void LBSReadAngularFluxes(std::shared_ptr<opensn::LBSProblem> lbs_problem,
                          const std::string& file_base);

} // namespace opensnlua
