// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/materials/interpolator/interpolator.h"
#include "framework/materials/interpolator/linear_interpolator.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/multi_group_xs/xsfile.h"
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#define XS_GETTER(method_name)                                                                     \
  py::cpp_function([](MultiGroupXS& self) { return convert_vector(self.method_name()); },          \
                   py::keep_alive<0, 1>())

namespace opensn
{

// Wrap multi-group cross section
void
WrapMultiGroupXS(py::module& xs)
{
  // clang-format off
  // sparse matrix
  auto sparse_matrix = py::class_<SparseMatrix>(
    xs,
    "SparseMatrix",
    R"(
    Sparse matrix wrapper.

    Wrapper of :cpp:class:`opensn::SparseMatrix`.
    )"
  );
  sparse_matrix.def(
    "GetValueIJ",
    &SparseMatrix::GetValueIJ,
    py::arg("i"),
    py::arg("j"),
    "Get the value at matrix entry ``(i, j)``."
  );
  sparse_matrix.def_property_readonly(
    "shape",
    [](const SparseMatrix& self) { return py::make_tuple(self.GetNumRows(), self.GetNumCols()); },
    "Get the matrix shape as ``(rows, columns)``."
  );
  sparse_matrix.def_property_readonly(
    "num_columns_per_row",
    [](const SparseMatrix& self)
    {
      std::vector<std::size_t> counts;
      counts.reserve(self.GetNumRows());
      for (const auto& column_indices : self.rowI_indices)
        counts.push_back(column_indices.size());
      return counts;
    },
    "Get the number of stored columns in each row."
  );

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
    This method mutates ``self`` by replacing its current contents.
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
    "Get total cross section."
  );
  multigroup_xs.def_property_readonly(
    "sigma_a",
    XS_GETTER(GetSigmaAbsorption),
    "Get absorption cross section."
  );
  multigroup_xs.def_property_readonly(
    "energy_deposition",
    XS_GETTER(GetEnergyDeposition),
    "Get energy deposition cross section."
  );
  multigroup_xs.def_property_readonly(
    "sigma_f",
    XS_GETTER(GetSigmaFission),
    "Get fission cross section."
  );
  multigroup_xs.def_property_readonly(
    "chi",
    XS_GETTER(GetChi),
    "Get neutron fission spectrum."
  );
  multigroup_xs.def_property_readonly(
    "nu_sigma_f",
    XS_GETTER(GetNuSigmaF),
    "Get neutron production due to fission."
  );
  multigroup_xs.def_property_readonly(
    "nu_prompt_sigma_f",
    XS_GETTER(GetNuPromptSigmaF),
    "Get prompt neutron production due to fission."
  );
  multigroup_xs.def_property_readonly(
    "nu_delayed_sigma_f",
    XS_GETTER(GetNuDelayedSigmaF),
    "Get delayed neutron production due to fission."
  );
  multigroup_xs.def_property_readonly(
    "production_matrix",
    [](const MultiGroupXS& self) { return self.GetProductionMatrix(); },
    "Get the fission production matrix."
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
    "Get inverse velocity."
  );
  multigroup_xs.def(
    "GetTransferMatrix",
    &MultiGroupXS::GetTransferMatrix,
    py::arg("ell"),
    py::return_value_policy::reference_internal,
    "Get the transfer matrix for moment ``ell``."
  );
  // clang-format on
}

// Wrap cross section interpolators
void
WrapInterpolator(py::module& xs)
{
  // clang-format off
  // xs type
  auto xs_type = py::enum_<XSType>(xs, "XSType");
  xs_type.value("Total",              XSType::Total);
  xs_type.value("Absorption",         XSType::Absorption);
  xs_type.value("Fission",            XSType::Fission);
  xs_type.value("NuFission",          XSType::NuFission);
  xs_type.value("Chi",                XSType::Chi);
  xs_type.value("ProductionMatrix",   XSType::ProductionMatrix);
  xs_type.value("NuPromptFission",    XSType::NuPromptFission);
  xs_type.value("NuDelayedFission",   XSType::NuDelayedFission);
  xs_type.value("Precursor",          XSType::Precursor);
  xs_type.value("Transfer",           XSType::Transfer);
  xs_type.value("Default",            XSType::Default);
  xs_type.def(
    "__or__",
    [](XSType a, XSType b) {
      return static_cast<XSType>(static_cast<std::uint64_t>(a) | static_cast<std::uint64_t>(b));
    }
  );

  // grid
  auto cartesian_grid = py::class_<CartesianGrid, std::shared_ptr<CartesianGrid>>(
    xs,
    "CartesianGrid",
    R"(
    Cartesian interpolation grid.

    Wrapper of :cpp:class:`opensn::CartesianGrid`.
    )"
  );
  cartesian_grid.def(
    py::init(
      [](const std::vector<std::vector<double>>& grid_data)
      {
        return std::make_shared<CartesianGrid>(grid_data);
      }),
    R"(
    Construct a Cartesian interpolation grid from a sequence of point arrays.

    Parameters
    ----------
    grid_data: Sequence[Sequence[float]]
        One sequence per dimension, each containing sorted unique grid points.
    )",
    py::arg("grid_data")
  );
  cartesian_grid.def(
    "__repr__",
    [](const CartesianGrid& self)
    {
      std::ostringstream os;
      os << "CartesianGrid([";

      const auto shape = self.GetShape();
      for (std::size_t d = 0; d < shape.size(); ++d)
      {
        if (d > 0)
          os << ", ";

        os << "[";
        const auto grid_d = self.GetGrid(static_cast<std::uint32_t>(d));
        for (std::size_t i = 0; i < grid_d.size(); ++i)
        {
          if (i > 0)
            os << ", ";
          os << grid_d[i];
        }
        os << "]";
      }

      os << "])";
      return os.str();
    }
  );

  // interpolator
  auto interpolator = py::class_<Interpolator, std::shared_ptr<Interpolator>>(
    xs,
    "Interpolator",
    R"(
    Cross-section interpolator.

    Wrapper of :cpp:class:`opensn::Interpolator`.
    )"
  );
  interpolator.def(
    "Evaluate",
    [](Interpolator& self, py::array_t<double, py::array::c_style | py::array::forcecast> state)
    {
      py::buffer_info buf = state.request();
      if (buf.ndim != 1)
        throw std::runtime_error("Interpolator.Evaluate expects a 1D NumPy array.");

      auto* ptr = static_cast<double*>(buf.ptr);
      const auto size = static_cast<std::size_t>(buf.shape[0]);
      return std::make_shared<MultiGroupXS>(self.Evaluate(std::span<double>(ptr, size)));
    },
    R"(
    Evaluate the interpolator at a state point.

    Parameters
    ----------
    state: numpy.ndarray
        One-dimensional NumPy array containing the interpolation state point.

    Returns
    -------
    pyopensn.xs.MultiGroupXS
        Interpolated multi-group cross section.
    )",
    py::arg("state")
  );

  // linear interpolator
  auto linear_interp = py::class_<LinearInterpolator, std::shared_ptr<LinearInterpolator>, Interpolator>(
    xs,
    "LinearInterpolator",
    R"(
    Cross-section linear interpolator.

    Wrapper of :cpp:class:`opensn::LinearInterpolator`.
    )"
  );
  linear_interp.def(
    py::init(
      [](const std::shared_ptr<CartesianGrid>& grid,
         py::array xs_data,
         XSType flag)
      {
        if (not grid)
          throw std::runtime_error("grid cannot be None.");

        const auto shape = grid->GetShape();
        if (xs_data.ndim() != static_cast<py::ssize_t>(shape.size()))
          throw std::runtime_error("xs_data rank does not match the CartesianGrid rank.");

        if (py::str(xs_data.dtype().attr("kind")).cast<std::string>() != "O")
          throw std::runtime_error("xs_data must be a NumPy object array of MultiGroupXS.");

        for (std::size_t d = 0; d < shape.size(); ++d)
        {
          if (xs_data.shape(static_cast<py::ssize_t>(d)) != static_cast<py::ssize_t>(shape[d]))
            throw std::runtime_error("xs_data shape does not match the CartesianGrid shape.");
        }

        py::buffer_info buf = xs_data.request();
        auto* base_ptr = static_cast<char*>(buf.ptr);
        std::vector<std::shared_ptr<MultiGroupXS>> xs_vec;
        xs_vec.reserve(grid->GetSize());
        std::vector<py::ssize_t> index(shape.size(), 0);

        for (std::size_t flat_i = 0; flat_i < grid->GetSize(); ++flat_i)
        {
          py::ssize_t byte_offset = 0;
          for (std::size_t d = 0; d < shape.size(); ++d)
            byte_offset += index[d] * buf.strides[d];

          // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
          auto* obj_ptr = reinterpret_cast<PyObject**>(base_ptr + byte_offset);
          xs_vec.push_back(
            py::reinterpret_borrow<py::object>(*obj_ptr).cast<std::shared_ptr<MultiGroupXS>>());

          for (std::size_t d = shape.size(); d-- > 0;)
          {
            ++index[d];
            if (index[d] < buf.shape[d])
              break;
            index[d] = 0;
          }
        }

        return std::make_shared<LinearInterpolator>(
          *grid, xs_vec, static_cast<std::uint64_t>(flag));
      }),
    R"(
    Construct a linear cross-section interpolator.

    Parameters
    ----------
    grid: pyopensn.xs.CartesianGrid
        Interpolation grid.
    xs_data: numpy.ndarray
        Object array of ``MultiGroupXS`` values whose shape matches ``grid``.
        The array is traversed in C-order before constructing the interpolator.
    flag: pyopensn.xs.XSType
        Bitmask describing which cross-section data are interpolated.
    )",
    py::arg("grid"),
    py::arg("xs_data"),
    py::arg("flag") = XSType::Default
  );

  // clang-format on
}

// Wrap the cross section components of OpenSn
void
py_xs(py::module& pyopensn)
{
  py::module xs = pyopensn.def_submodule("xs", "Cross section module.");
  WrapMultiGroupXS(xs);
  WrapInterpolator(xs);
}

} // namespace opensn
