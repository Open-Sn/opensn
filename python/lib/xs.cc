// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/material_property.h"
#include <memory>
#include <string>
#include <vector>

#define XS_GETTER(method_name) [](MultiGroupXS& self) { return convert_vector(self.method_name()); }

namespace opensn
{

// clang-format off

// Wrap multi-group cross section
void wrap_multigroup_xs(py::module &xs)
{
  // multi-group cross section
  auto multigroup_xs = py::class_<MultiGroupXS, std::shared_ptr<MultiGroupXS>, MaterialProperty>(
    xs, "MultiGroupXS",
    R"(
    Multi-group cross section.

    Wrapper of :cpp:class:`opensn::MultiGroupXS`.
    )");

  multigroup_xs.def(
    py::init(
      [](void)
      {
        std::shared_ptr<MultiGroupXS> xs = std::make_shared<MultiGroupXS>();
        multigroup_xs_stack.push_back(xs);
        return xs;
      }),
    "Create an empty multi-group cross section.");

  multigroup_xs.def(
    "Initialize",
    [](MultiGroupXS &self, double sigma_t, double c) { self.Initialize(sigma_t, c); },
    "Makes a simple material with a 1-group cross-section set.");

  multigroup_xs.def("SetScalingFactor", &MultiGroupXS::SetScalingFactor, "Scale the cross sections by the specified factor.");
  multigroup_xs.def_property_readonly("num_groups", &MultiGroupXS::GetNumGroups, "Get number of energy groups.");
  multigroup_xs.def_property_readonly("scattering_order", &MultiGroupXS::GetScatteringOrder, "Get Legendre scattering order.");
  multigroup_xs.def_property_readonly("num_precursors", &MultiGroupXS::GetNumPrecursors, "Get number of precursors.");
  multigroup_xs.def_property_readonly("is_fissionable", &MultiGroupXS::IsFissionable, "Check if the material is fissile.");
  multigroup_xs.def_property_readonly("scaling_factor", &MultiGroupXS::GetScalingFactor, "Get the arbitrary scaling factor.");

  multigroup_xs.def_property_readonly("sigma_t", XS_GETTER(GetSigmaTotal), "Get total cross section.", py::keep_alive<0, 1>());
  multigroup_xs.def_property_readonly("sigma_a", XS_GETTER(GetSigmaAbsorption), "Get absorption cross section.", py::keep_alive<0, 1>());
  multigroup_xs.def_property_readonly("sigma_f", XS_GETTER(GetSigmaFission), "Get fission cross section.", py::keep_alive<0, 1>());
  multigroup_xs.def_property_readonly("chi", XS_GETTER(GetChi), "Get neutron fission spectrum.", py::keep_alive<0, 1>());
  multigroup_xs.def_property_readonly("nu_sigma_f", XS_GETTER(GetNuSigmaF), "Get neutron production due to fission.", py::keep_alive<0, 1>());
  multigroup_xs.def_property_readonly("nu_prompt_sigma_f", XS_GETTER(GetNuPromptSigmaF), "Get prompt neutron production due to fission.", py::keep_alive<0, 1>());
  multigroup_xs.def_property_readonly("nu_delayed_sigma_f", XS_GETTER(GetNuDelayedSigmaF), "Get delayed neutron production due to fission.", py::keep_alive<0, 1>());
  multigroup_xs.def_property_readonly("inv_velocity", XS_GETTER(GetInverseVelocity), "Get inverse velocity.", py::keep_alive<0, 1>());
}

// Wrap create and load cross sections
void wrap_create_load(py::module &xs)
{
  xs.def(
    "CreateSimpleOneGroup",
    [](double sigma_t, double c)
    {
      std::shared_ptr<MultiGroupXS> xs = std::make_shared<MultiGroupXS>();
      xs->Initialize(sigma_t, c);
      return xs;
    },
    "Create a one-group cross section.");

  xs.def(
    "LoadFromOpenSn",
    [](const std::string &file_name)
    {
      std::shared_ptr<MultiGroupXS> xs = std::make_shared<MultiGroupXS>();
      xs->Initialize(file_name);
      return xs;
    },
    "Load multi-group cross sections from an OpenSn cross section input file.");

  xs.def(
    "LoadFromOpenMC",
    [](const std::string &file_name, const std::string &dataset_name, double temperature)
    {
      std::shared_ptr<MultiGroupXS> xs = std::make_shared<MultiGroupXS>();
      xs->Initialize(file_name, dataset_name, temperature);
      return xs;
    },
    "Load multi-group cross sections from an OpenMC cross-section file.");
}

// Wrap the cross section components of OpenSn
void py_xs(py::module &pyopensn)
{
  py::module xs = pyopensn.def_submodule("xs", "Cross section module.");
  wrap_multigroup_xs(xs);
  wrap_create_load(xs);
}

// clang-format on

} // namespace opensn
