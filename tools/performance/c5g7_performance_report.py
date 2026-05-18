#!/usr/bin/env python3
"""Run C5G7 k-eigenvalue solver variants and summarize convergence/runtime metrics."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EXE = REPO_ROOT / "build" / "python" / "opensn"
TRANSPORT_KEIGEN_DIR = (
    REPO_ROOT / "test" / "python" / "modules" / "linear_boltzmann_solvers" / "transport_keigen"
)

DEFAULT_METHODS = ("pi", "pi_cmfd", "pi_scdsa", "pi_smm", "jfnk")


@dataclass
class RunMetrics:
    method: str
    returncode: int
    elapsed_seconds: float | None = None
    k_eff: float | None = None
    sweeps: int | None = None
    func_evals: int | None = None
    cmfd_metrics: dict[str, list[float]] | None = None
    output: str = ""


def parse_elapsed_seconds(output: str) -> float | None:
    match = re.search(r"Elapsed execution time:\s+(\d+):(\d+):(\d+(?:\.\d+)?)", output)
    if not match:
        return None
    hours, minutes, seconds = match.groups()
    return int(hours) * 3600.0 + int(minutes) * 60.0 + float(seconds)


def parse_k_eff(output: str) -> float | None:
    for pattern in (
        r"Python k-eigenvalue:\s*([-+0-9.eE]+)",
        r"(?:PI|NLKE) final,\s+status\s+=\s+\w+,\s+k_eff\s+=\s*([-+0-9.eE]+)",
    ):
        match = re.search(pattern, output)
        if match:
            return float(match.group(1))
    return None


def parse_sweeps(output: str) -> int | None:
    match = re.search(r"Python sweeps:\s*(\d+)", output)
    if match:
        return int(match.group(1))
    match = re.search(r"sweeps\s*=\s*(\d+)", output)
    return int(match.group(1)) if match else None


def parse_func_evals(output: str) -> int | None:
    match = re.search(r"func evals\s*=\s*(\d+)", output)
    return int(match.group(1)) if match else None


def parse_cmfd_metrics(output: str) -> dict[str, list[float]]:
    metrics: dict[str, list[float]] = {}
    for line in output.splitlines():
        if "CMFD_ACCEL" not in line:
            continue
        fields = {}
        for token in line.split():
            if "=" in token:
                key, value = token.split("=", 1)
                fields[key] = value
        category = fields.get("category")
        metric = fields.get("metric")
        value = fields.get("value")
        if category is not None and metric is not None and value is not None:
            metrics.setdefault(f"{category}.{metric}", []).append(float(value))
    return metrics


def parse_run(method: str, returncode: int, output: str) -> RunMetrics:
    return RunMetrics(
        method=method,
        returncode=returncode,
        elapsed_seconds=parse_elapsed_seconds(output),
        k_eff=parse_k_eff(output),
        sweeps=parse_sweeps(output),
        func_evals=parse_func_evals(output),
        cmfd_metrics=parse_cmfd_metrics(output),
        output=output,
    )


def py_assignment(name: str, value: str) -> str:
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'{name}="{escaped}"'


def build_command(
    exe: Path,
    ranks: int,
    mesh_type: str,
    method: str,
) -> list[str]:
    cmd = [
        "mpirun",
        "-np",
        str(ranks),
        str(exe),
        "-i",
        "c5g7/c5g7.py",
        "--suppress-color",
        "--py",
        py_assignment("mesh_type", mesh_type),
        "--py",
        py_assignment("k_method", method),
    ]
    if method == "pi_cmfd":
        cmd.extend(
            [
                "--py",
                "cmfd_verbose=True",
                "--py",
                "cmfd_group_aggregation_size=4",
            ]
        )
    return cmd


def run_method(
    exe: Path,
    ranks: int,
    mesh_type: str,
    method: str,
) -> RunMetrics:
    cmd = build_command(exe, ranks, mesh_type, method)
    print(f"Running c5g7:{mesh_type}:{method}: {' '.join(cmd)}", flush=True)
    completed = subprocess.run(
        cmd,
        cwd=TRANSPORT_KEIGEN_DIR,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return parse_run(method, completed.returncode, completed.stdout)


def fmt_float(value: float | None, precision: int = 6) -> str:
    return "n/a" if value is None else f"{value:.{precision}g}"


def fmt_seconds(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1f}s"


def work_count(metric: RunMetrics) -> str:
    if metric.sweeps is not None:
        return str(metric.sweeps)
    if metric.func_evals is not None:
        return str(metric.func_evals)
    return "n/a"


def work_label(metric: RunMetrics) -> str:
    if metric.sweeps is not None:
        return "sweeps"
    if metric.func_evals is not None:
        return "evals"
    return "n/a"


def relative_time(reference: RunMetrics | None, metric: RunMetrics) -> float | None:
    if (
        reference is None
        or reference.elapsed_seconds is None
        or metric.elapsed_seconds in (None, 0)
    ):
        return None
    return reference.elapsed_seconds / metric.elapsed_seconds


def k_difference(reference: RunMetrics | None, metric: RunMetrics) -> float | None:
    if reference is None or reference.k_eff is None or metric.k_eff is None:
        return None
    return abs(reference.k_eff - metric.k_eff)


def metric_value(metrics: dict[str, list[float]] | None, key: str, default: str = "n/a") -> str:
    if metrics is None:
        return default
    values = metrics.get(key, [])
    if not values:
        return default
    value = (
        max(values)
        if key.endswith((".skipped", ".invalid_k", ".nonfinite_flux"))
        else values[-1]
    )
    return fmt_float(value)


def print_report(metrics: list[RunMetrics]) -> None:
    baseline = next((metric for metric in metrics if metric.method == "pi"), None)
    cmfd = next((metric for metric in metrics if metric.method == "pi_cmfd"), None)
    header = (
        "method",
        "k_eff",
        "k_diff_vs_pi",
        "work",
        "work_type",
        "time",
        "speedup_vs_pi",
        "speedup_vs_cmfd",
        "cmfd_unknowns",
        "cmfd_skipped",
    )
    print("\nC5G7 performance report")
    print(" ".join(f"{column:>16}" for column in header))
    for metric in metrics:
        cmfd_speedup = (
            metric.elapsed_seconds / cmfd.elapsed_seconds
            if (
                cmfd is not None
                and metric.elapsed_seconds is not None
                and cmfd.elapsed_seconds not in (None, 0)
            )
            else None
        )
        values = (
            metric.method,
            fmt_float(metric.k_eff, 8),
            fmt_float(k_difference(baseline, metric), 3),
            work_count(metric),
            work_label(metric),
            fmt_seconds(metric.elapsed_seconds),
            fmt_float(relative_time(baseline, metric), 3),
            fmt_float(cmfd_speedup, 3),
            metric_value(metric.cmfd_metrics, "coarse_mesh.global_unknowns"),
            metric_value(metric.cmfd_metrics, "correction.skipped"),
        )
        print(" ".join(f"{value:>16}" for value in values))


def write_outputs(output_dir: Path | None, metrics: list[RunMetrics]) -> None:
    if output_dir is None:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    for metric in metrics:
        output_path = output_dir / f"c5g7_{metric.method}.out"
        output_path.write_text(metric.output, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--ranks", type=int, default=32)
    parser.add_argument("--mesh-type", choices=("coarse", "fine"), default="coarse")
    parser.add_argument("--method", choices=DEFAULT_METHODS, action="append")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    selected_methods = args.method or list(DEFAULT_METHODS)
    metrics = [
        run_method(
            args.exe,
            args.ranks,
            args.mesh_type,
            method,
        )
        for method in selected_methods
    ]
    write_outputs(args.output_dir, metrics)

    failed = [metric for metric in metrics if metric.returncode != 0]
    if failed:
        print(
            "C5G7 report failed for: " + ", ".join(metric.method for metric in failed),
            file=sys.stderr,
        )
        return 1

    print_report(metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
