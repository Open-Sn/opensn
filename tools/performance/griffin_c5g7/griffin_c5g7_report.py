#!/usr/bin/env python3
"""Run c5g7_griffin_sweep.py over the Griffin direction counts and compare per-sweep times.

Griffin reference values are for one NVIDIA A100 (40 GB), first-order DFEM, 4 threads per
element and 64 per team (N. Choi et al., Nucl. Sci. Eng. 2026, Fig. 6a and Sec. IV.A). Only the
384-direction time (0.169 s) is stated in the text; the others are read from the log-scale
Fig. 6a and are uncertain by roughly 10%. Griffin grind time uses the same normalization as
OpenSn's sweep_time_per_unknown: sweep time / (elements * 4 * 7 groups * directions).
"""

from __future__ import annotations

import argparse
import re
import statistics
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
INPUT = SCRIPT_DIR / "c5g7_griffin_sweep.py"

NUM_ELEMENTS = 101728
NUM_GROUPS = 7
DOFS_PER_ELEMENT = 4
GRIFFIN_A100_SWEEP_TIME = {24: 0.030, 48: 0.040, 96: 0.057, 192: 0.095, 384: 0.169}

TIMING_RE = re.compile(
    r"avg_sweep_time = ([-+0-9.eE]+) s, sweep_time_per_unknown = ([-+0-9.eE]+) ns")
UNKNOWNS_RE = re.compile(r"unknowns = (\d+), lagged_unknowns = (\d+)")


def griffin_grind_ns(num_dirs: int) -> float:
    unknowns = NUM_ELEMENTS * DOFS_PER_ELEMENT * NUM_GROUPS * num_dirs
    return GRIFFIN_A100_SWEEP_TIME[num_dirs] * 1.0e9 / unknowns


def py_assignment(name: str, value) -> str:
    if isinstance(value, str):
        escaped = value.replace("\\", "\\\\").replace('"', '\\"')
        return f'{name}="{escaped}"'
    return f"{name}={value}"


def run_case(args, num_dirs: int, sweep_type: str, log_dir: Path) -> dict:
    params = {
        "num_dirs": num_dirs,
        "n_polar": args.n_polar,
        "sweep_type": sweep_type,
        "angle_agg": args.angle_agg,
        "use_gpus": not args.cpu,
        "sweeps_per_outer": args.sweeps_per_outer,
        "num_outers": args.num_outers,
        "mesh_file": str(args.mesh_file.resolve()),
        "vacuum_boundaries": args.vacuum,
    }
    cmd = args.launcher.split() + ["-n", str(args.ranks)] if args.launcher else []
    cmd += [str(args.exe), "-i", str(INPUT), "--suppress-color"]
    for name, value in params.items():
        cmd += ["--py", py_assignment(name, value)]
    if args.extra:
        cmd += ["--py", f"extra_problem_params={args.extra}"]
    if args.extra_options:
        cmd += ["--py", f"extra_options={args.extra_options}"]

    tag = f"{args.label}_{sweep_type}_{num_dirs}"
    print(f"[{tag}] {' '.join(cmd)}", flush=True)
    completed = subprocess.run(cmd, cwd=log_dir, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, check=False)
    (log_dir / f"{tag}.log").write_text(completed.stdout)

    timings = [(float(a), float(b)) for a, b in TIMING_RE.findall(completed.stdout)]
    unknowns = UNKNOWNS_RE.findall(completed.stdout)
    result = {"tag": tag, "num_dirs": num_dirs, "sweep_type": sweep_type,
              "returncode": completed.returncode, "num_solves": len(timings)}
    # The first power iteration includes one-time GPU setup, so it is never used. A case without
    # at least one more completed power iteration has no steady-state time.
    steady = timings[1:]
    if steady:
        result["sweep_time"] = statistics.median(t for t, _ in steady)
        result["grind_ns"] = statistics.median(g for _, g in steady)
    if unknowns:
        result["unknowns"], result["lagged"] = (int(v) for v in unknowns[-1])
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exe", type=Path, default=REPO_ROOT / "build" / "python" / "opensn")
    parser.add_argument("--mesh-file", type=Path, default=SCRIPT_DIR / "c5g7_griffin.msh")
    parser.add_argument("--launcher", default="mpiexec",
                        help="MPI launcher; empty string to run the executable directly")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--dirs", type=int, nargs="+", default=sorted(GRIFFIN_A100_SWEEP_TIME))
    parser.add_argument("--sweep-types", nargs="+", default=["AAH", "CBC"])
    parser.add_argument("--n-polar", type=int, default=4)
    parser.add_argument("--angle-agg", default="polar")
    parser.add_argument("--sweeps-per-outer", type=int, default=10)
    parser.add_argument("--num-outers", type=int, default=3)
    parser.add_argument("--extra-options", default="",
                        help="Python dict literal merged into the problem's options block, "
                             "e.g. \"{'fused_aao_sweep': True}\"")
    parser.add_argument("--extra", default="",
                        help="Python dict literal merged into the problem parameters")
    parser.add_argument("--cpu", action="store_true", help="run without use_gpus")
    parser.add_argument("--vacuum", action="store_true",
                        help="diagnostic: vacuum instead of reflecting xmin/ymin boundaries")
    parser.add_argument("--label", default="opensn")
    parser.add_argument("--log-dir", type=Path, default=Path.cwd() / "griffin_c5g7_logs")
    args = parser.parse_args()

    if args.num_outers < 2:
        parser.error("--num-outers must be at least 2 because the first power iteration is a "
                     "warm-up")
    if not args.mesh_file.is_file():
        print(f"{args.mesh_file} not found; run make_c5g7_quad_mesh.py first", file=sys.stderr)
        return 1
    # Cases run with the log directory as their working directory.
    args.exe = args.exe.resolve()
    args.log_dir.mkdir(parents=True, exist_ok=True)

    # Alternate the sweep-type order from one direction count to the next, so drift during the
    # run does not favor one sweep type. Cases run, and are logged, in this order.
    cases = []
    for i, d in enumerate(args.dirs):
        sweep_types = args.sweep_types if i % 2 == 0 else list(reversed(args.sweep_types))
        cases += [(d, st) for st in sweep_types]
    results = [run_case(args, d, st, args.log_dir) for d, st in cases]

    print()
    print(f"{'case':<24} {'dirs':>5} {'sweep (s)':>11} {'grind (ns)':>11} "
          f"{'Griffin (s)':>12} {'Griffin (ns)':>13} {'ratio':>7} {'lagged':>8}")
    status = 0
    for r in results:
        if r["returncode"] != 0 or "sweep_time" not in r:
            print(f"{r['tag']:<24} {r['num_dirs']:>5}  FAILED (rc={r['returncode']}, "
                  f"{r['num_solves']} timed solves, see {args.log_dir / (r['tag'] + '.log')})")
            status = 1
            continue
        d = r["num_dirs"]
        expected = NUM_ELEMENTS * DOFS_PER_ELEMENT * NUM_GROUPS * d
        if r.get("unknowns") != expected:
            # A different problem size cannot be compared with Griffin's fixed reference.
            print(f"{r['tag']:<24} {d:>5}  INVALID ({r.get('unknowns')} unknowns, Griffin has "
                  f"{expected})")
            status = 1
            continue
        g_time = GRIFFIN_A100_SWEEP_TIME.get(d)
        g_grind = griffin_grind_ns(d) if g_time else None
        ratio = r["sweep_time"] / g_time if g_time else None
        print(f"{r['tag']:<24} {d:>5} {r['sweep_time']:>11.4e} {r['grind_ns']:>11.4f} "
              f"{g_time if g_time else float('nan'):>12.3f} "
              f"{g_grind if g_grind else float('nan'):>13.4f} "
              f"{ratio if ratio else float('nan'):>7.2f} {r.get('lagged', 0):>8}")
    print("\nratio = OpenSn sweep time / Griffin A100 sweep time (< 1 means OpenSn is faster).")
    return status


if __name__ == "__main__":
    sys.exit(main())
