// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"

#ifdef OPENSN_WITH_PRECICE

#include "adapters/precice/cell_power.h"
#include "adapters/precice/precice_adapter.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/utils/collective_error.h"
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <optional>
#include <vector>

namespace opensn
{

void
WrapPreciceAdapter(py::module& coupling)
{
  // clang-format off
  auto adapter = py::class_<coupling::PreciceAdapter>(
    coupling,
    "PreciceAdapter",
    R"(
    One-way OpenSn preCICE volume-coupling participant.

    The adapter registers locally owned OpenSn cell centroids and sends
    cell-integrated power in watts. It contains no ghost cells and performs no
    file-based transfer.
    )"
  );
  adapter.def(
    py::init<std::filesystem::path, std::string, std::string, std::string, double>(),
    R"(
    Construct an OpenSn preCICE participant.

    Collective over OpenSn's communicator, which is passed to preCICE.

    Parameters
    ----------
    configuration_file : pathlib.Path or str
        preCICE XML configuration file.
    participant_name : str, default='OpenSn'
        Participant name declared in the preCICE configuration.
    mesh_name : str, default='OpenSn-Cells'
        Coupling mesh provided by OpenSn. It must be 3D for 1D and 3D OpenSn meshes and may be
        2D or 3D for 2D OpenSn meshes.
    power_data_name : str, default='TransportCellPower'
        Scalar cell-integrated power data, in watts per cell.
    length_scale : float, default=1.0
        Factor converting OpenSn mesh coordinates to the receiver's length unit, for example
        ``0.01`` for an OpenSn mesh in centimeters and a receiver in meters. Must be finite and
        positive. Transferred power is unchanged.
    )",
    py::arg("configuration_file"),
    py::arg("participant_name") = "OpenSn",
    py::arg("mesh_name") = "OpenSn-Cells",
    py::arg("power_data_name") = "TransportCellPower",
    py::arg("length_scale") = 1.0
  );
  adapter.def(
    "WriteSteadyCellPower",
    [](coupling::PreciceAdapter& self,
       const FieldFunctionGridBased& field_function,
       const double source_rate,
       const double joules_per_energy_unit,
       const std::size_t component,
       const std::optional<double> transverse_measure,
       const std::vector<unsigned int>& block_ids,
       const py::object& averaging_group)
    {
      auto cell_data = coupling::MakeCellPowerData(field_function,
                                                   component,
                                                   source_rate,
                                                   joules_per_energy_unit,
                                                   transverse_measure,
                                                   block_ids);
      if (not averaging_group.is_none())
      {
        std::vector<std::int64_t> group_ids;
        group_ids.reserve(cell_data.GetCells().size());
        RunWithCollectiveErrors(
          [&]
          {
            for (const auto& cell : cell_data.GetCells())
              group_ids.push_back(
                averaging_group(cell.block_id,
                                cell.centroid.x,
                                cell.centroid.y,
                                cell.centroid.z).cast<std::int64_t>());
          });
        cell_data.ConservativelyAverageByGroup(group_ids);
      }
      py::gil_scoped_release release;
      self.WriteSteadyCellPower(cell_data);
    },
    R"(
    Integrate a deposition field and send fixed power until coupling ends.

    The field is the energy deposited per unit volume, with volume in the OpenSn mesh length unit,
    for the source of the transport input, for example from
    ``problem.CreateFieldFunction(name, "energy_deposition")``. Each owned-cell integral is
    multiplied by ``source_rate``, ``joules_per_energy_unit``, and, for 1D and 2D meshes,
    ``transverse_measure`` to produce watts per cell. The mesh must be Cartesian. The method may be
    called once per adapter. Collective. See the preCICE coupling section of the user guide.

    Parameters
    ----------
    field_function : pyopensn.fieldfunc.FieldFunctionGridBased
        Grid field containing volumetric energy deposition.
    source_rate : float
        Physical source particles per second divided by the input's total source rate for the
        physical body. For 1D and 2D meshes, that input rate is the rate per unit area or depth,
        such as ``production_rate + inflow_rate`` from ``ComputeBalanceTable()``, times
        ``transverse_measure``. Must be finite and nonnegative.
    joules_per_energy_unit : float, default=1.602176634e-13
        Joules per field energy unit. The default converts MeV to joules.
    component : int, default=0
        Field component to integrate.
    transverse_measure : float, optional
        Required for 1D and 2D meshes and rejected for 3D meshes. It is the cross-sectional area
        of a 1D mesh or the out-of-plane depth of a 2D mesh, in the OpenSn mesh length unit.
    block_ids : List[int], default=[]
        Blocks whose cells are coupled. Empty couples every cell. Each listed block must contain at
        least one cell. Power in other blocks is not transferred.
    averaging_group : Callable[[int, float, float, float], int], optional
        Map ``(block_id, x, y, z)``, with the centroid in OpenSn mesh coordinates, to a global
        control-volume ID. Within each ID, power density is replaced by its global volume-weighted
        average. This conservative projection can remove finite-element undershoots. Use a modest
        number of IDs, not one per cell.
    )",
    py::arg("field_function"),
    py::arg("source_rate"),
    py::arg("joules_per_energy_unit") = 1.602176634e-13,
    py::arg("component") = 0,
    py::arg("transverse_measure") = py::none(),
    py::arg("block_ids") = std::vector<unsigned int>{},
    py::arg("averaging_group") = py::none()
  );
  // clang-format on
}

void
py_coupling(py::module& pyopensn)
{
  py::module coupling = pyopensn.def_submodule("coupling", "Multiphysics coupling module.");
  WrapPreciceAdapter(coupling);
}

} // namespace opensn

#endif // OPENSN_WITH_PRECICE
