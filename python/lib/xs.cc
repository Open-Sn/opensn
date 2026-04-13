// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/multi_group_xs/xsfile.h"
#include <pybind11/stl.h>
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

    The Python API currently has two types of methods:

    - Creation/loading methods such as ``CreateSimpleOneGroup``,
      ``LoadFromOpenSn``, and ``LoadFromOpenMC`` populate an existing object.
    - ``Scale`` mutates the current object.
    - ``Combine`` returns a new cross-section object and does not mutate inputs.
    )"
  );
  multigroup_xs.def(
    py::init(
      []()
      {
        return std::make_shared<MultiGroupXS>();
      }),
    "Create an empty multi-group cross section."
  );
  multigroup_xs.def(
    "CreateSimpleOneGroup",
    [](MultiGroupXS& self, double sigma_t, double c, double velocity) {
      self = MultiGroupXS::CreateSimpleOneGroup(sigma_t, c, velocity);
    },
    R"(
    Populate this object with a one-group cross section.

    Parameters
    ----------
    sigma_t: float
        Total cross section.
    c: float
        Scattering ratio.
    velocity: float, optional
        Group velocity. If provided and positive, inverse velocity
        is populated with 1.0/velocity.

    Notes
    -----
    This method mutates ``self`` by replacing its current contents.
    )",
    py::arg("sigma_t"),
    py::arg("c"),
    py::arg("velocity") = 0.0
  );
  multigroup_xs.def(
    "LoadFromOpenSn",
    [](MultiGroupXS& self, const std::string& file_name)
    {
      self = MultiGroupXS::LoadFromOpenSn(file_name);
    },
    py::arg("file_name"),
    R"(
    Load multi-group cross sections from an OpenSn cross section input file
    into this object.

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
       M_GFROM_GTO_VAL 0 0 0 value
       .
       M_GFROM_GTO_VAL moment gfrom gto value
       .
       M_GFROM_GTO_VAL nmom-1 ng-1 ng-1 value
       TRANSFER_MOMENTS_END

    Notes
    -----
    This method mutates ``self`` by replacing its current contents.
    )"
  );
  multigroup_xs.def_static(
    "Combine",
    &MultiGroupXS::Combine,
    R"(
    Return a new combined cross-section.

    Parameters
    ----------

    combinations: List[Tuple[pyopensn.xs.MultiGroupXS, float]]
        List of ``(cross_section, density)`` pairs.
        The density values are linear weights used to combine raw cross sections.

    Returns
    -------
    pyopensn.xs.MultiGroupXS
        A new combined cross-section object. The input cross sections are not
        modified.

    Notes
    -----
    Let :math:`d_i` be the supplied density for cross section :math:`i`.

    - Raw XS terms are density-weighted sums:
      :math:`\sigma = \sum_i d_i \sigma_i`
      (e.g. total, absorption, fission, transfer, production).
    - Named custom 1D XS are preserved and combined with the same density weighting.
    - Fission spectra and precursor yields are weighted by fissile density
      fraction so their sums remain normalized.
    - All inputs must have the same number of groups.
    - If inverse velocity is present, all inputs must have identical values.

    Examples
    --------

    >>> xs_1 = MultiGroupXS()
    >>> xs_1.CreateSimpleOneGroup(sigma_t=1, c=0.5)
    >>> xs_2 = MultiGroupXS()
    >>> xs_2.CreateSimpleOneGroup(sigma_t=2, c=1./3.)
    >>> combo = [
    ...     ( xs_1, 0.5 ),
    ...     ( xs_2, 3.0 )
    ... ]
    >>> xs_combined = MultiGroupXS.Combine(combo)
    )",
    py::arg("combinations")
  );
  multigroup_xs.def(
    "LoadFromOpenMC",
    [](MultiGroupXS& self,
       const std::string& file_name,
       const std::string& dataset_name,
       double temperature,
       const std::vector<std::string>& extra_xs_names)
    {
      self = MultiGroupXS::LoadFromOpenMC(file_name, dataset_name, temperature, extra_xs_names);
    },
    R"(
    Load multi-group cross sections from an OpenMC cross-section file into this
    object.

    Notes
    -----
    When delayed-neutron data is present in the OpenMC MGXS library, this
    method also imports:

    - the number of precursor groups
    - prompt and delayed fission production data
    - precursor decay constants
    - delayed emission spectra

    If delayed data is present but ``chi-prompt`` is not available, the steady
    fission spectrum ``chi`` is used as the prompt spectrum.
    )",
    py::arg("file_name"),
    py::arg("dataset_name"),
    py::arg("temperature"),
    py::arg("extra_xs_names") = std::vector<std::string>()
  );
  multigroup_xs.def(
    "LoadFromCEPXS",
    [](MultiGroupXS& self, const std::string& file_name, int material_id)
    {
      self = MultiGroupXS::LoadFromCEPXS(file_name, material_id);
    },
    "Load multi-group cross sections from a CEPXS cross-section file.",
    py::arg("file_name"),
    py::arg("material_id") = 0
  );
  multigroup_xs.def(
    "Scale",
    &MultiGroupXS::Scale,
    R"(
    Scale the cross sections in-place.

    Notes
    -----
    Scaling does not compound. Each call scales from the original baseline data.
    Named custom 1D XS are scaled along with the standard 1D cross-section data.
    )",
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
  multigroup_xs.def(
    "GetScaleFactor",
    &MultiGroupXS::GetScaleFactor,
    "Get the scaling factor."
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
    "energy_deposition",
    XS_GETTER(GetEnergyDeposition),
    "Get energy deposition cross section.",
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
  multigroup_xs.def(
    "has_custom_xs",
    &MultiGroupXS::HasCustomXS,
    "Check if a custom XS is available.",
    py::arg("name")
  );
  multigroup_xs.def(
    "get_custom_xs",
    [](MultiGroupXS& self, const std::string& name)
    { return convert_vector(self.GetCustomXS(name)); },
    "Get a custom XS vector.",
    py::arg("name")
  );
  multigroup_xs.def(
    "custom_xs_names",
    &MultiGroupXS::GetCustomXSNames,
    "Get a list of custom XS entries."
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
