// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "python/lib/console.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <pybind11/pybind11.h>
#include <memory>

using namespace opensn;
namespace py = pybind11;

namespace unit_tests
{

void data_types_Test00();
void TestCFunction();
void math_Test01_WDD_IJK_Sweep();
void math_Test02_ParallelVector();
void math_SDM_Test01_Continuous(std::shared_ptr<MeshContinuum> grid,
                                std::string sdm_type,
                                bool export_vtk);
void math_SDM_Test02_Discontinuous(std::shared_ptr<MeshContinuum> grid,
                                   std::string sdm_type,
                                   bool export_vtk);
void SimTest01_FV(std::shared_ptr<MeshContinuum> grid);
void SimTest02_FV(std::shared_ptr<MeshContinuum> grid);
void SimTest03_PWLC(std::shared_ptr<MeshContinuum> grid);
void SimTest04_PWLC(std::shared_ptr<MeshContinuum> grid);
void SimTest06_WDD(std::shared_ptr<MeshContinuum> grid);
void SimTest91_PWLD(std::shared_ptr<MeshContinuum> grid);
void SimTest93_RayTracing(std::shared_ptr<MeshContinuum> grid);
void acceleration_Diffusion_CFEM(std::shared_ptr<MeshContinuum> grid);
void acceleration_Diffusion_DFEM(std::shared_ptr<MeshContinuum> grid);

} // namespace unit_tests
