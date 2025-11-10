// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "caribou/caribou.h"

namespace crb = caribou;

namespace opensn
{

// Wrap device settings
void
WrapDeviceSettings(py::module& device)
{
  // clang-format off
  device.def(
    "get_device_count",
    []() {
      return crb::get_num_gpus();
    },
    "Get the number of GPU devices visible to the current MPI rank."
  );
  device.def(
    "get_current_device",
    []() {
      return crb::get_current_device();
    },
    "Get the current GPU device number for the current MPI rank."
  );
  device.def(
    "set_device",
    [](int device_num) {
      crb::set_device(device_num);
    },
    R"(
    Set GPU device for the current MPI rank using device number

    The device number must range from 0 up to (but not including) ``get_device_count()``.
    )",
    py::arg("device_num")
  );
  // clang-format on
}

// Wrap the device setting components of OpenSn
void
py_device(py::module& pyopensn)
{
  py::module device = pyopensn.def_submodule("device", "GPU settings module.");
  WrapDeviceSettings(device);
}

} // namespace opensn
