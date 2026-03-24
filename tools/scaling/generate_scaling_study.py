#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from pathlib import Path
from lib import generate_strong_scaling, generate_weak_scaling

base_dir = Path(__file__).resolve().parent

geo_filename = base_dir / "lib/cube.geo"
"""Name of the Gmsh geometry file."""

gmsh_binary = "gmsh"
"""Path to the Gmsh binary. Default is 'gmsh' assuming it's in the system PATH."""

opensn_binary = Path(__file__).resolve().parents[2] / "build/python/opensn"
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
"""List of node counts to generate SLURM job files for."""

ncores = 96
"""Number of CPU cores per node."""

partition = ""
"""Name of the SLURM partition to submit jobs to."""

ngpus = 0
"""Number of GPUs per node (default: 0 for CPU-only)."""

gpu_config = "--gpus-per-task=1"
"""GPU configuration for SLURM job (default: --gpus-per-task=1)."""

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
    args = parser.parse_args()

    inputs = {
        "opensn_binary": opensn_binary,
        "gmsh_binary": gmsh_binary,
        "geo_filename": geo_filename,
        "ncores": ncores,
        "partition": partition,
        "ngpus": ngpus,
        "gpu_config": gpu_config
    }
    if args.type == "strong":
        generate_strong_scaling(nodes, divisor=strong_divisor, **inputs)
    else:
        generate_weak_scaling(nodes, weak_divisors, **inputs)
