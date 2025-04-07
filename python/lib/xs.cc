// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include <memory>
#include <string>
#include <vector>

#define XS_GETTER(method_name) [](MultiGroupXS& self) { return convert_vector(self.method_name()); }

namespace opensn
{

// Wrap multi-group cross section
void
WrapMultiGroupXS(py::module& xs)
{
  // clang-format off
  // multi-group cross section
  auto multigroup_xs = py::class_<MultiGroupXS, std::shared_ptr<MultiGroupXS>>(
    xs,
    "MultiGroupXS",
    R"(
    Multi-group cross section.

    Wrapper of :cpp:class:`opensn::MultiGroupXS`.
    )"
  );
  multigroup_xs.def(
    py::init(
      []()
      {
        std::shared_ptr<MultiGroupXS> xs = std::make_shared<MultiGroupXS>();
        multigroup_xs_stack.push_back(xs);
        return xs;
      }),
    "Create an empty multi-group cross section."
  );
  multigroup_xs.def(
    "CreateSimpleOneGroup",
    [](MultiGroupXS& self, double sigma_t, double c) {
      self.Initialize(sigma_t, c);
    },
    R"(
    Create a one-group cross section.

    Parameters
    ----------
    sigma_t: float
        Total cross section.
    c: float
        Scattering ratio.
    )",
    py::arg("sigma_t"),
    py::arg("c")
  );
  multigroup_xs.def(
    "LoadFromOpenSn",
    [](MultiGroupXS& self, const std::string& file_name)
    {
      self.Initialize(file_name);
    },
    py::arg("file_name"),
    R"(
    Load multi-group cross sections from an OpenSn cross section input file.

    Format is as follows (for transfers, gprime denotes the departing group and g is the arrival
    group).

    .. code-block:: none

       # Add comment lines, as needed
       NUM_GROUPS ng
       NUM_MOMENTS nmom

       SIGMA_T_BEGIN
       0 value
       .
       .
       ng-1 value
       SIGMA_T_END

       SIGMA_A_BEGIN
       0 value
       .
       .
       ng-1 value
       SIGMA_A_END

       TRANSFER_MOMENTS_BEGIN
       M_GPRIME_G_VAL 0 0 0 value
       .
       M_GPRIME_G_VAL moment gprime g value
       .
       M_GPRIME_G_VAL nmom-1 ng-1 ng-1 value
       TRANSFER_MOMENTS_END
    )"
  );
  multigroup_xs.def(
    "LoadFromOpenMC",
    [](MultiGroupXS& self, const std::string& file_name, const std::string& dataset_name,
       double temperature)
    {
      self.Initialize(file_name, dataset_name, temperature);
    },
    "Load multi-group cross sections from an OpenMC cross-section file.",
    py::arg("file_name"),
    py::arg("dataset_name"),
    py::arg("temperature")
  );
  multigroup_xs.def(
    "SetScalingFactor",
    &MultiGroupXS::SetScalingFactor,
    "Scale the cross sections by the specified factor.",
    py::arg("factor")
  );
  multigroup_xs.def_property_readonly(
    "num_groups",
    &MultiGroupXS::GetNumGroups,
    "Get number of energy groups."
  );
  multigroup_xs.def_property_readonly(
    "scattering_order",
    &MultiGroupXS::GetScatteringOrder,
    "Get Legendre scattering order."
  );
  multigroup_xs.def_property_readonly(
    "num_precursors",
    &MultiGroupXS::GetNumPrecursors,
    "Get number of precursors."
  );
  multigroup_xs.def_property_readonly(
    "is_fissionable",
    &MultiGroupXS::IsFissionable,
    "Check if the material is fissile."
  );
  multigroup_xs.def_property_readonly(
    "scaling_factor",
    &MultiGroupXS::GetScalingFactor,
    "Get the arbitrary scaling factor."
  );
  multigroup_xs.def_property_readonly(
    "sigma_t",
    XS_GETTER(GetSigmaTotal),
    "Get total cross section.",
    py::keep_alive<0, 1>()
  );
  multigroup_xs.def_property_readonly(
    "sigma_a",
    XS_GETTER(GetSigmaAbsorption),
    "Get absorption cross section.",
    py::keep_alive<0, 1>()
  );
  multigroup_xs.def_property_readonly(
    "sigma_f",
    XS_GETTER(GetSigmaFission),
    "Get fission cross section.",
    py::keep_alive<0, 1>()
  );
  multigroup_xs.def_property_readonly(
    "chi",
    XS_GETTER(GetChi),
    "Get neutron fission spectrum.",
    py::keep_alive<0, 1>()
  );
  multigroup_xs.def_property_readonly(
    "nu_sigma_f",
    XS_GETTER(GetNuSigmaF),
    "Get neutron production due to fission.",
    py::keep_alive<0, 1>()
  );
  multigroup_xs.def_property_readonly(
    "nu_prompt_sigma_f",
    XS_GETTER(GetNuPromptSigmaF),
    "Get prompt neutron production due to fission.",
    py::keep_alive<0, 1>()
  );
  multigroup_xs.def_property_readonly(
    "nu_delayed_sigma_f",
    XS_GETTER(GetNuDelayedSigmaF),
    "Get delayed neutron production due to fission.",
    py::keep_alive<0, 1>()
  );
  multigroup_xs.def_property_readonly(
    "inv_velocity",
    XS_GETTER(GetInverseVelocity),
    "Get inverse velocity.",
    py::keep_alive<0, 1>()
  );
  // clang-format on
}

// Wrap the cross section components of OpenSn
void
py_xs(py::module& pyopensn)
{
  py::module xs = pyopensn.def_submodule("xs", "Cross section module.");
  WrapMultiGroupXS(xs);
}

} // namespace opensn
