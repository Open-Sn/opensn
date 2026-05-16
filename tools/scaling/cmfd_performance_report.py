#!/usr/bin/env python3
"""Run PI, accelerated PI, and NLKE k-eigenvalue benchmarks and summarize metrics."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EXE = REPO_ROOT / "build" / "python" / "opensn"
TRANSPORT_KEIGEN_DIR = (
    REPO_ROOT / "test" / "python" / "modules" / "linear_boltzmann_solvers" / "transport_keigen"
)


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    ranks: int
    baseline: str
    cmfd: str
    nlke: str | None = None
    scdsa: str | None = None
    smm: str | None = None
    cmfd_args: tuple[str, ...] = ("cmfd_verbose=True",)
    nlke_args: tuple[str, ...] = ("solver_name='nlke'",)
    scdsa_args: tuple[str, ...] = ("solver_name='scdsa'",)
    smm_args: tuple[str, ...] = ("solver_name='smm'",)


@dataclass
class RunMetrics:
    name: str
    label: str
    returncode: int
    elapsed_seconds: float | None = None
    k_eff: float | None = None
    sweeps: int | None = None
    func_evals: int | None = None
    cmfd_metrics: dict[str, list[float]] = field(default_factory=dict)
    output: str = ""


CASES = {
    "unstructured8": BenchmarkCase(
        name="unstructured8",
        ranks=8,
        baseline="keigenvalue_transport_3d_2g_unstructured_8rank.py",
        cmfd="keigenvalue_transport_3d_2g_unstructured_8rank_cmfd_agg.py",
        nlke="keigenvalue_transport_3d_2g_unstructured_8rank.py",
        scdsa="keigenvalue_transport_3d_2g_unstructured_8rank.py",
        smm="keigenvalue_transport_3d_2g_unstructured_8rank.py",
    ),
}


def parse_elapsed_seconds(output: str) -> float | None:
    match = re.search(r"Elapsed execution time:\s+(\d+):(\d+):(\d+(?:\.\d+)?)", output)
    if not match:
        return None
    hours, minutes, seconds = match.groups()
    return int(hours) * 3600.0 + int(minutes) * 60.0 + float(seconds)


def parse_k_eff(output: str) -> float | None:
    match = re.search(r"Python k-eigenvalue:\s*([-+0-9.eE]+)", output)
    return float(match.group(1)) if match else None


def parse_sweeps(output: str) -> int | None:
    match = re.search(r"Python sweeps:\s*(\d+)", output)
    if match:
        return int(match.group(1))

    final = re.search(r"sweeps\s*=\s*(\d+)", output)
    return int(final.group(1)) if final else None


def parse_func_evals(output: str) -> int | None:
    match = re.search(r"func evals\s*=\s*(\d+)", output)
    return int(match.group(1)) if match else None


def parse_cmfd_metrics(output: str) -> dict[str, list[float]]:
    metrics: dict[str, list[float]] = {}
    for line in output.splitlines():
        if "CMFD_METRIC" not in line:
            continue
        fields = {}
        for token in line.split():
            if "=" in token:
                key, value = token.split("=", 1)
                fields[key] = value
        category = fields.get("c")
        metric = fields.get("m")
        value = fields.get("v")
        if category is None or metric is None or value is None:
            continue
        metrics.setdefault(f"{category}.{metric}", []).append(float(value))
    return metrics


def parse_run(name: str, label: str, returncode: int, output: str) -> RunMetrics:
    return RunMetrics(
        name=name,
        label=label,
        returncode=returncode,
        elapsed_seconds=parse_elapsed_seconds(output),
        k_eff=parse_k_eff(output),
        sweeps=parse_sweeps(output),
        func_evals=parse_func_evals(output),
        cmfd_metrics=parse_cmfd_metrics(output),
        output=output,
    )


def build_command(exe: Path, ranks: int, input_file: str, py_args: tuple[str, ...]) -> list[str]:
    cmd = ["mpirun", "-np", str(ranks), str(exe), "-i", input_file, "--suppress-color"]
    for arg in py_args:
        cmd.extend(["--py", arg])
    return cmd


def run_case(exe: Path, case: BenchmarkCase, label: str, input_file: str, py_args: tuple[str, ...]):
    cmd = build_command(exe, case.ranks, input_file, py_args)
    print(f"Running {case.name}:{label}: {' '.join(cmd)}", flush=True)
    completed = subprocess.run(
        cmd,
        cwd=TRANSPORT_KEIGEN_DIR,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return parse_run(case.name, label, completed.returncode, completed.stdout)


def fmt_float(value: float | None, precision: int = 6) -> str:
    return "n/a" if value is None else f"{value:.{precision}g}"


def fmt_seconds(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1f}s"


def metric_value(metrics: dict[str, list[float]], key: str, default: str = "n/a") -> str:
    values = metrics.get(key, [])
    if not values:
        return default
    value = max(values) if key.endswith((".skipped", ".invalid_k", ".nonfinite_flux")) else values[-1]
    return fmt_float(value)


def method_work(metric: RunMetrics | None) -> str:
    if metric is None:
        return "n/a"
    if metric.sweeps is not None:
        return str(metric.sweeps)
    if metric.func_evals is not None:
        return str(metric.func_evals)
    return "n/a"


def method_time(metric: RunMetrics | None) -> str:
    return fmt_seconds(metric.elapsed_seconds if metric is not None else None)


def relative_time(reference: RunMetrics, metric: RunMetrics | None) -> float | None:
    if metric is None or reference.elapsed_seconds is None or metric.elapsed_seconds in (None, 0):
        return None
    return reference.elapsed_seconds / metric.elapsed_seconds


def k_difference(reference: RunMetrics, metric: RunMetrics | None) -> float | None:
    if metric is None or reference.k_eff is None or metric.k_eff is None:
        return None
    return abs(reference.k_eff - metric.k_eff)


def print_report(rows: list[tuple[RunMetrics, RunMetrics, RunMetrics | None, RunMetrics | None, RunMetrics | None]]):
    header = (
        "case",
        "pi_sweeps",
        "cmfd_sweeps",
        "nlke_evals",
        "scdsa_sweeps",
        "smm_p0_sweeps",
        "sweep_speedup",
        "pi_time",
        "cmfd_time",
        "nlke_time",
        "scdsa_time",
        "smm_p0_time",
        "time_speedup",
        "cmfd_vs_nlke",
        "cmfd_vs_scdsa",
        "cmfd_vs_smm_p0",
        "k_diff",
        "nlke_k_diff",
        "scdsa_k_diff",
        "smm_p0_k_diff",
        "coarse_unknowns",
        "undersized",
        "avg_faces",
        "max_faces",
        "skipped",
    )
    print("\nCMFD performance report")
    print(" ".join(f"{column:>16}" for column in header))
    for baseline, cmfd, nlke, scdsa, smm in rows:
        sweep_speedup = (
            baseline.sweeps / cmfd.sweeps
            if baseline.sweeps is not None and cmfd.sweeps not in (None, 0)
            else None
        )
        k_diff = (
            abs(baseline.k_eff - cmfd.k_eff)
            if baseline.k_eff is not None and cmfd.k_eff is not None
            else None
        )
        values = (
            baseline.name,
            method_work(baseline),
            method_work(cmfd),
            method_work(nlke),
            method_work(scdsa),
            method_work(smm),
            fmt_float(sweep_speedup, 3),
            method_time(baseline),
            method_time(cmfd),
            method_time(nlke),
            method_time(scdsa),
            method_time(smm),
            fmt_float(relative_time(baseline, cmfd), 3),
            fmt_float(relative_time(nlke, cmfd) if nlke is not None else None, 3),
            fmt_float(relative_time(scdsa, cmfd) if scdsa is not None else None, 3),
            fmt_float(relative_time(smm, cmfd) if smm is not None else None, 3),
            fmt_float(k_diff, 3),
            fmt_float(k_difference(cmfd, nlke), 3),
            fmt_float(k_difference(cmfd, scdsa), 3),
            fmt_float(k_difference(cmfd, smm), 3),
            metric_value(cmfd.cmfd_metrics, "coarse_mesh.global_unknowns"),
            metric_value(cmfd.cmfd_metrics, "coarse_mesh.undersized_coarse_cells"),
            metric_value(cmfd.cmfd_metrics, "coarse_mesh.average_faces_per_coarse_cell"),
            metric_value(cmfd.cmfd_metrics, "coarse_mesh.max_faces_per_coarse_cell"),
            metric_value(cmfd.cmfd_metrics, "correction.skipped"),
        )
        print(" ".join(f"{value:>16}" for value in values))


def write_outputs(output_dir: Path | None, metrics: list[RunMetrics]):
    if output_dir is None:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    for metric in metrics:
        output_path = output_dir / f"{metric.name}_{metric.label}.out"
        output_path.write_text(metric.output, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--case", choices=sorted(CASES), action="append")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    selected_cases = [CASES[name] for name in (args.case or ["unstructured8"])]
    all_metrics: list[RunMetrics] = []
    report_rows: list[tuple[RunMetrics, RunMetrics, RunMetrics | None, RunMetrics | None, RunMetrics | None]] = []
    for case in selected_cases:
        baseline = run_case(args.exe, case, "pi", case.baseline, ())
        cmfd = run_case(args.exe, case, "cmfd", case.cmfd, case.cmfd_args)
        nlke = (
            run_case(args.exe, case, "nlke", case.nlke, case.nlke_args)
            if case.nlke is not None
            else None
        )
        scdsa = (
            run_case(args.exe, case, "scdsa", case.scdsa, case.scdsa_args)
            if case.scdsa is not None
            else None
        )
        smm = (
            run_case(args.exe, case, "smm_p0", case.smm, case.smm_args)
            if case.smm is not None
            else None
        )
        all_metrics.extend(metric for metric in (baseline, cmfd, nlke, scdsa, smm) if metric is not None)
        report_rows.append((baseline, cmfd, nlke, scdsa, smm))
        failed = any(
            metric is not None and metric.returncode != 0
            for metric in (baseline, cmfd, nlke, scdsa, smm)
        )
        if failed:
            write_outputs(args.output_dir, all_metrics)
            print(f"Case {case.name} failed; rerun with --output-dir to retain output.", file=sys.stderr)
            return 1

    print_report(report_rows)
    write_outputs(args.output_dir, all_metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
