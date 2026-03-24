#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate GMSH input files for scaling test with unstructured mesh.
"""

from pathlib import Path
from jinja2 import Environment, FileSystemLoader
import os
import subprocess
import shutil

base_dir = Path(__file__).resolve().parent
env = Environment(
    loader=FileSystemLoader(base_dir),
    autoescape=False
)
slurm_template = env.get_template("slurm_template.txt")
py_template = env.get_template("unstructured.py")


def run_gmsh(gmsh_binary, input_geo, divisor, output_msh):
    """Run Gmsh on a .geo file to generate a mesh."""

    cmd = [
        gmsh_binary,
        "-3",
        "-v", "0",
        "-setnumber", "divisor", str(divisor),
        "-o", str(output_msh),
        str(input_geo),
    ]
    print("Generating mesh with command: {}".format(" ".join(cmd)))
    subprocess.run(cmd, check=True)


def make_slurm_keys(base_name, opensn_binary, n_cores, partition, ngpus, gpu_config):
    """Create a dictionary of keys for SLURM job file generation."""
    keys = {
        "base_name": base_name,
        "opensn_binary": opensn_binary,
        "partition": partition,
    }
    if ngpus > 0:
        keys["use_gpus"] = True
        keys["gpu_options"] = f"#SBATCH {gpu_config}"
        keys["n_tasks"] = ngpus
        keys["n_cores"] = n_cores // ngpus
    else:
        keys["use_gpus"] = False
        keys["gpu_options"] = ""
        keys["n_tasks"] = n_cores
        keys["n_cores"] = 1
    return keys


def create_slurm_file(keys, n_nodes, input_script, prefix, output_dir):
    """Create a SLURM job file to run OpenSn."""

    slurm_content = slurm_template.render(**keys, n_nodes=n_nodes, input_script=input_script)
    slurm_filename = output_dir / f"{prefix}_{n_nodes}.sh"
    print(f"Generating SLURM file {slurm_filename}")
    with open(slurm_filename, "w") as slurm_file:
        slurm_file.write(slurm_content)
    return slurm_filename


def generate_strong_scaling(nodes, **kwargs):
    """Generate files for strong scaling study."""

    # copy necessary files to output directory
    out_dir = base_dir.parent / "output" / "strong"
    os.makedirs(out_dir, exist_ok=True)
    shutil.copyfile(
        base_dir / "xs_168g.xs",
        out_dir / "xs_168g.xs"
    )

    # generate the mesh
    run_gmsh(
        kwargs["gmsh_binary"],
        kwargs["geo_filename"],
        kwargs["divisor"],
        output_msh=out_dir / "strong_scaling.msh"
    )

    # create SLURM job files
    keys = make_slurm_keys(
        base_name="strong",
        opensn_binary=kwargs["opensn_binary"],
        n_cores=kwargs["ncores"],
        partition=kwargs["partition"],
        ngpus=kwargs["ngpus"],
        gpu_config=kwargs["gpu_config"]
    )
    script_name = out_dir / "strong_scaling.py"
    slurm_scripts = []
    for n_nodes in nodes:
        slurm_script = create_slurm_file(
            keys=keys,
            n_nodes=n_nodes,
            input_script=script_name,
            prefix="strong",
            output_dir=out_dir
        )
        slurm_scripts.append(slurm_script)
    with open(out_dir / "submit_jobs.sh", "w") as launch_file:
        launch_file.write("#!/bin/bash\n\n")
        for slurm_script in slurm_scripts:
            launch_file.write(f"sbatch {slurm_script.name}\n")

    # generate the Python script
    print("Generating strong scaling Python script.")
    script_content = py_template.render(**keys, mesh_file="strong_scaling.msh")
    with open(script_name, "w") as script_file:
        script_file.write(script_content)


def generate_weak_scaling(nodes, divisors, **kwargs):
    """Generate files for weak scaling study."""

    # copy necessary files to output directory
    out_dir = base_dir.parent / "output" / "weak"
    os.makedirs(out_dir, exist_ok=True)
    shutil.copyfile(
        base_dir / "xs_168g.xs",
        out_dir / "xs_168g.xs"
    )

    # generate the meshes
    for n_nodes, divisor in zip(nodes, divisors):
        run_gmsh(
            kwargs["gmsh_binary"],
            kwargs["geo_filename"],
            divisor,
            output_msh=out_dir / f"weak_scaling_{n_nodes}.msh"
        )

    # create SLURM job files
    keys = make_slurm_keys(
        base_name="weak",
        opensn_binary=kwargs["opensn_binary"],
        n_cores=kwargs["ncores"],
        partition=kwargs["partition"],
        ngpus=kwargs["ngpus"],
        gpu_config=kwargs["gpu_config"]
    )
    slurm_scripts = []
    for n_nodes in nodes:
        slurm_script = create_slurm_file(
            keys=keys,
            n_nodes=n_nodes,
            input_script=out_dir / f"weak_scaling_{n_nodes}.py",
            prefix="weak",
            output_dir=out_dir
        )
        slurm_scripts.append(slurm_script)
    with open(out_dir / "submit_jobs.sh", "w") as launch_file:
        launch_file.write("#!/bin/bash\n\n")
        for slurm_script in slurm_scripts:
            launch_file.write(f"sbatch {slurm_script.name}\n")

    # generate the Python script
    print("Generating weak scaling Python script.")
    for n_nodes in nodes:
        script_content = py_template.render(**keys, mesh_file=f"weak_scaling_{n_nodes}.msh")
        script_name = out_dir / f"weak_scaling_{n_nodes}.py"
        with open(script_name, "w") as script_file:
            script_file.write(script_content)
