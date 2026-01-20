// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_app.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include "test/python/src/bindings.h"
#include "petsc.h"

namespace mpi = mpicpp_lite;

void
register_bindings()
{
  // clang-format off
  opensnpy::Console::BindModule([](py::module& m) {
    m.def("data_types_Test00", &unit_tests::data_types_Test00); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("TestCFunction", &unit_tests::TestCFunction); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("math_Test01_WDD_IJK_Sweep", &unit_tests::math_Test01_WDD_IJK_Sweep); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("math_Test02_ParallelVector", &unit_tests::math_Test02_ParallelVector); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("math_SDM_Test01_Continuous", &unit_tests::math_SDM_Test01_Continuous); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("math_SDM_Test02_Discontinuous", &unit_tests::math_SDM_Test02_Discontinuous); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("SimTest01_FV", &unit_tests::SimTest01_FV); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("SimTest02_FV", &unit_tests::SimTest02_FV); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("SimTest03_PWLC", &unit_tests::SimTest03_PWLC); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("SimTest04_PWLC", &unit_tests::SimTest04_PWLC); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("SimTest06_WDD", &unit_tests::SimTest06_WDD); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("SimTest91_PWLD", &unit_tests::SimTest91_PWLD); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("SimTest93_RayTracing", &unit_tests::SimTest93_RayTracing); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("acceleration_Diffusion_CFEM", &unit_tests::acceleration_Diffusion_CFEM); });

  opensnpy::Console::BindModule([](py::module& m) {
    m.def("acceleration_Diffusion_DFEM", &unit_tests::acceleration_Diffusion_DFEM); });
  // clang-format on
}

int
main(int argc, char** argv)
{
  mpi::Environment env(argc, argv, mpi::ThreadSupport::MULTIPLE);

  PetscCall(PetscInitializeNoArguments());

  int retval = EXIT_SUCCESS;
  try
  {
    py::scoped_interpreter guard{};
    register_bindings();
    opensnpy::PyApp app(MPI_COMM_WORLD);
    retval = app.Run(argc, argv);
  }
  catch (...)
  {
    std::fprintf(stderr, "Unknown fatal error\n");
    retval = EXIT_FAILURE;
  }

  PetscFinalize();

  return retval;
}
