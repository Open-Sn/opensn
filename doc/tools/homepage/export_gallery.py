# SPDX-License-Identifier: MIT
"""Export the C5G7 thermal flux used by the README artwork.

Run from the repository root with the OpenSn console::

    build/python/opensn -i doc/tools/homepage/export_gallery.py

Override gallery_output to choose the destination (default: gallery-output).
"""

from pathlib import Path
import runpy

if "opensn_console" not in globals():
    raise RuntimeError("Run this input with the OpenSn console, not Python.")

root = Path(__file__).resolve().parents[3]
api = globals()
output = Path(globals().get("gallery_output", "gallery-output")).resolve()
if api["rank"] == 0:
    output.mkdir(parents=True, exist_ok=True)
api["MPIBarrier"]()

source = root / "test/python/modules/linear_boltzmann_solvers/transport_keigen/c5g7/c5g7.py"
state = runpy.run_path(
    str(source),
    init_globals={**globals(), "mesh_type": "coarse", "k_method": "pi_scdsa"},
    run_name="__main__",
)
solver = state["k_solver"]
if not solver.IsConverged():
    raise RuntimeError("The C5G7 solve did not converge; do not use its output as artwork.")
# Same discretization and acceptance as the existing C5G7 regression.
if abs(solver.GetEigenvalue() - 1.192559) > 1.0e-4:
    raise RuntimeError("C5G7 does not match the existing regression reference.")
api["FieldFunctionGridBased"].ExportMultipleToPVTU(
    [state["phys"].GetScalarFluxFieldFunction()[6]], str(output / "c5g7_thermal")
)
