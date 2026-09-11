#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import glob
import statistics
import warnings
import yaml
import matplotlib.pyplot as plt
from datetime import datetime
from argparse import ArgumentParser
from matplotlib.ticker import NullLocator
from pathlib import Path
from generate_scaling_study import extra_data


def extract_data(filename):
    """Extract n, average sweep time, and its standard deviation from a file.

    Each outer repetition in the output file reports its own
    ``avg_sweep_time``. The first repetition is always discarded (warm-up),
    and the mean/sample standard deviation of the remaining sweep times are
    used as the metric and its error.
    """

    match = re.search(r"_(\d+)\.out$", filename)
    if not match:
        return None
    n = int(match.group(1))

    sweep_times = []

    avg_time_re = re.compile(r"avg_sweep_time\s*=\s*([0-9.eE+-]+)\s*s")

    with open(filename, "r") as f:
        for line in f:
            avg_match = avg_time_re.search(line)
            if avg_match:
                sweep_times.append(float(avg_match.group(1)))

    # discard the first (warm-up) iteration
    sweep_times = sweep_times[1:]

    if not sweep_times:
        return None

    avg_time = statistics.mean(sweep_times)
    std_time = statistics.stdev(sweep_times) if len(sweep_times) > 1 else 0.0

    return n, avg_time, std_time


def compute_efficiency(data):
    """Compute efficiency and its error relative to the baseline (first) run.

    Efficiency is ``base_time * 100 / t``, a ratio of two independent
    measurements. Its error is obtained via standard error propagation for a
    quotient: ``(err/eff)^2 = (base_std/base_time)^2 + (std/t)^2``. The
    baseline point has no independent second measurement to divide by, so
    only its own measurement uncertainty contributes: ``err = eff *
    (base_std/base_time)``.
    """

    base_time, base_std = data[0][1], data[0][2]
    base_rel_var = (base_std / base_time) ** 2
    efficiency = []
    error = []
    for i, (_, t, std) in enumerate(data):
        eff = base_time * 100.0 / t
        if i == 0:
            err = eff * (base_std / base_time)
        else:
            err = eff * (base_rel_var + (std / t) ** 2) ** 0.5
        efficiency.append(eff)
        error.append(err)
    return efficiency, error


def plot_data(data, output_file, with_history):
    """Plot the data and save to a file."""

    n_nodes = [d[0] for d in data]
    efficiency, error = compute_efficiency(data)

    history = {}
    if with_history and (Path(__file__).resolve().parent / "history.yaml").exists():
        with open("history.yaml", "r") as f:
            history_dict = yaml.safe_load(f)
        name = extra_data["name"]
        history_label = f"{name}_weak_scaling"
        if history_dict is not None and history_label in history_dict:
            history_data = history_dict[history_label]
            history["nodes"] = history_data["nodes"]
            history["efficiency"] = history_data["efficiency"]
            history["error"] = history_data.get("error", [0.0] * len(history["nodes"]))

    fig, ax = plt.subplots()
    ax.errorbar(n_nodes, efficiency, yerr=error, marker="o", color="xkcd:cerulean",
                label="efficiency", capsize=3)
    xticks = n_nodes.copy()
    if history:
        ax.errorbar(history["nodes"], history["efficiency"], yerr=history["error"],
                    marker="o", color="xkcd:coral", label="history", capsize=3)
        xticks = sorted(set(n_nodes) | set(history["nodes"]))
    elif with_history:
        warnings.warn(
            "History file not found or history label not in file. "
            "Plotting without history."
        )
    ax.set_xlabel("Number of nodes")
    ax.set_xscale("log")
    ax.set_xticks(xticks, xticks)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(bottom=0.0, top=max(e + err for e, err in zip(efficiency, error)) + 10.0)
    ax.set_ylabel("Efficiency (%)")
    ax.set_title("Node-to-node weak scaling")
    ax.grid(True, which="both")
    ax.legend()
    fig.savefig(output_file)
    plt.show()


def export_data(data, output_file):
    """Export data to a YAML file."""

    efficiency, error = compute_efficiency(data)

    name = extra_data["name"]
    label = f"{name}_weak_scaling"
    export_dict = None
    if (Path(__file__).resolve().parent / "history.yaml").exists():
        with open(output_file, "r") as f:
            export_dict = yaml.safe_load(f)
    if export_dict is None:
        export_dict = {}
    export_dict[label] = {
        "description": extra_data["description"],
        "time": datetime.now().isoformat(),
        "nodes": [d[0] for d in data],
        "efficiency": efficiency,
        "error": error
    }
    with open(output_file, "w") as f:
        yaml.dump(export_dict, f)


if __name__ == "__main__":

    # read command-line arguments
    parser = ArgumentParser(description="Plot scaling data from output files.")
    parser.add_argument(
        "--output",
        type=str,
        default="weak_scaling_plot.pdf",
        help="Filename for the output plot (default: weak_scaling_plot.pdf)."
    )
    parser.add_argument(
        "--dir",
        type=str,
        default="output/weak_cpu",
        help="Folder to find weak scaling result (default: output/weak_cpu)."
    )
    parser.add_argument(
        "--history",
        type=str,
        choices=["none", "comp", "save"],
        default="none",
        help=(
            "History mode for the plot: "
            "none (only plot current data), "
            "comp (compare with history in the same plot without saving), "
            "or save (plot and overwrite current history value). "
            "(default: none)"
        ),
    )
    args = parser.parse_args()

    # get files matching the prefix in the input directory
    input_dir = Path(__file__).resolve().parent / args.dir
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory {input_dir} does not exist.")
    files = glob.glob(f"{input_dir}/weak_*.out")
    if not files:
        raise FileNotFoundError(f"No files found matching weak_*.out in {input_dir}")

    # extract sweep time
    data = []
    for f in files:
        result = extract_data(f)
        if result:
            data.append(result)
    if not data:
        raise ValueError("No valid data found.")
    data.sort(key=lambda x: x[0])

    # plot
    with_history = (args.history == "comp")
    plot_data(data, args.output, with_history)

    # export data to YAML
    if args.history == "save":
        export_data(data, "history.yaml")
