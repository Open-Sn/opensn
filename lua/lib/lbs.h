// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"

namespace opensnlua
{

std::vector<std::shared_ptr<opensn::FieldFunctionGridBased>>
LBSGetScalarFieldFunctionList(std::shared_ptr<opensn::LBSSolver> lbs_solver);

std::map<std::string, std::vector<double>>
LBSComputeLeakage(std::shared_ptr<opensn::DiscreteOrdinatesSolver> solver,
                  const std::vector<std::string>& bnd_names);

void LBSWriteFluxMoments(std::shared_ptr<opensn::LBSSolver> lbs_solver,
                         const std::string& file_base);

void LBSCreateAndWriteSourceMoments(std::shared_ptr<opensn::LBSSolver> lbs_solver,
                                    const std::string& file_base);

void LBSReadFluxMomentsAndMakeSourceMoments(std::shared_ptr<opensn::LBSSolver> lbs_solver,
                                            const std::string& file_base,
                                            bool single_file_flag);

void LBSReadSourceMoments(std::shared_ptr<opensn::LBSSolver> lbs_solver,
                          const std::string& file_base,
                          bool single_file_flag);

void LBSReadFluxMoments(std::shared_ptr<opensn::LBSSolver> lbs_solver,
                        const std::string& file_base,
                        bool single_file_flag);

void LBSWriteAngularFluxes(std::shared_ptr<opensn::LBSSolver> lbs_solver,
                           const std::string& file_base);

void LBSReadAngularFluxes(std::shared_ptr<opensn::LBSSolver> lbs_solver,
                          const std::string& file_base);

} // namespace opensnlua
