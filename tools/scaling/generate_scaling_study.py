#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from pathlib import Path
from lib import generate_strong_scaling, generate_weak_scaling

base_dir = Path(__file__).resolve().parent

environment = ""
"""Command to activate environement script."""

geo_filename = base_dir / "lib/cube.geo"
"""Name of the Gmsh geometry file."""

gmsh_binary = "gmsh"
"""Path to the Gmsh binary. Default is 'gmsh' assuming it's in the system PATH."""

opensn_binary = base_dir.parents[1] / "build/python/opensn"
"""Path to the OpenSn binary. Default is '../../../build/python/opensn'."""

strong_divisor = 39
"""Divisor for Gmsh to control mesh resolution for strong scaling study (default: 39)."""

# 64 tasks per node with 256 cells per task
weak_divisors = [15, 19, 24, 31, 39, 49, 62, 78, 98, 124]
# 96 tasks per node with 256 cell per task
# weak_divisors = [17, 22, 28, 35, 44, 56, 71, 89, 113, 137]
# 64 tasks per node with 2048 cells per task
# weak_divisors = [30, 39, 48, 62, 78, 98, 124]

nodes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
"""List of node counts to generate job files for."""

ncores = 96
"""
Number of CPU cores available per node.

If this value is not evenly divisible by the number of GPUs, round it down.
"""

partition = ""
"""Name of the partition/queue to submit jobs to."""

ngpus = 0
"""Number of GPUs per node (default: 0 for CPU-only)."""

gpu_config = {
    "slurm": "#SBATCH --gpus-per-task=1",
    "lc-flux": "#flux: --setattr=gpumode=SPX",
    "alcf-pbs": ""
}
"""GPU configurations."""

user_defined_config = ""
"""
Custom user-defined argument for launching with GPU.

If this variable is defined, it will overwrite default options in ``gpu_config``.
"""

extra_data = {
    "name": "label_name",
    "description": "Description of the cluster and scaling study."
}

if __name__ == "__main__":

    parser = ArgumentParser(
        description="Generate files for strong or weak scaling studies with OpenSn.")
    parser.add_argument(
        "--type",
        type=str,
        choices=["strong", "weak"],
        help="Type of scaling test to generate files for."
    )
    parser.add_argument(
        "--processor",
        type=str,
        choices=["cpu", "gpu"],
        default="cpu",
        help="Processor target for the generated study files. Defaults to CPU."
    )
    parser.add_argument(
        "--engine",
        type=str,
        choices=["slurm", "lc-flux", "alcf-pbs"],
        default="slurm",
        help="Job submitting system. Defaults to slurm."
    )
    args = parser.parse_args()

    if args.processor == "gpu" and ngpus == 0:
        raise ValueError("Please specify the number of GPUs per node.")

    if user_defined_config:
        gpu_option = user_defined_config
    else:
        gpu_option = gpu_config[args.engine]

    inputs = {
        "opensn_binary": opensn_binary,
        "gmsh_binary": gmsh_binary,
        "geo_filename": geo_filename,
        "ncores": ncores,
        "partition": partition,
        "processor": args.processor,
        "ngpus": ngpus,
        "engine": args.engine,
        "environment": environment,
        "gpu_config": gpu_option
    }
    if args.type == "strong":
        generate_strong_scaling(nodes, divisor=strong_divisor, **inputs)
    else:
        generate_weak_scaling(nodes, weak_divisors, **inputs)
