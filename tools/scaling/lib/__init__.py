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
templates = {
    "slurm": env.get_template("slurm_template.txt"),
    "lc-flux": env.get_template("flux_template.txt"),
    "alcf-pbs": env.get_template("pbs_template.txt")
}
commands = {
    "slurm": "sbatch",
    "lc-flux": "flux batch",
    "alcf-pbs": "qsub"
}
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


def make_launch_script_keys(
    base_name,
    outdir,
    opensn_binary,
    processor,
    n_cores,
    partition,
    ngpus,
    environment,
    gpu_config
):
    """Create a dictionary of keys for job file generation."""
    keys = {
        "base_name": base_name,
        "outdir": str(outdir),
        "opensn_binary": opensn_binary,
        "environment": environment,
        "partition": partition,
        "processor": processor
    }
    if processor == "gpu":
        keys["use_gpus"] = True
        keys["gpu_options"] = gpu_config
        keys["n_tasks"] = ngpus
        keys["n_cores"] = n_cores // ngpus
    else:
        keys["use_gpus"] = False
        keys["gpu_options"] = ""
        keys["n_tasks"] = n_cores
        keys["n_cores"] = 1
    return keys


def create_job_script(engine, keys, n_nodes, input_script, prefix, output_dir):
    """Create a job file to run OpenSn."""

    content = templates[engine].render(**keys, n_nodes=n_nodes, input_script=input_script)
    fname = output_dir / f"{prefix}_{n_nodes}.sh"
    print(f"Generating job file {fname}")
    with open(fname, "w") as job_file:
        job_file.write(content)
    return fname


def generate_strong_scaling(nodes, **kwargs):
    """Generate files for strong scaling study."""

    # copy necessary files to output directory
    processor = kwargs["processor"]
    engine = kwargs["engine"]
    out_dir = base_dir.parent / "output" / f"strong_{processor}"
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

    # create job files
    keys = make_launch_script_keys(
        base_name=f"strong_{processor}",
        outdir=out_dir,
        opensn_binary=kwargs["opensn_binary"],
        processor=processor,
        n_cores=kwargs["ncores"],
        partition=kwargs["partition"],
        ngpus=kwargs["ngpus"],
        environment=kwargs["environment"],
        gpu_config=kwargs["gpu_config"]
    )
    script_name = out_dir / "strong_scaling.py"
    scripts = []
    for n_nodes in nodes:
        script = create_job_script(
            engine=engine,
            keys=keys,
            n_nodes=n_nodes,
            input_script=script_name,
            prefix="strong",
            output_dir=out_dir
        )
        scripts.append(script)
    with open(out_dir / "submit_jobs.sh", "w") as launch_file:
        launch_file.write("#!/bin/bash\n\n")
        for script in scripts:
            launch_file.write(f"{commands[engine]} {script.name}\n")

    # generate the Python script
    print("Generating strong scaling Python script.")
    script_content = py_template.render(**keys, mesh_file="strong_scaling.msh")
    with open(script_name, "w") as script_file:
        script_file.write(script_content)


def generate_weak_scaling(nodes, divisors, **kwargs):
    """Generate files for weak scaling study."""

    # copy necessary files to output directory
    processor = kwargs["processor"]
    engine = kwargs["engine"]
    out_dir = base_dir.parent / "output" / f"weak_{processor}"
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

    # create job files
    keys = make_launch_script_keys(
        base_name=f"weak_{processor}",
        outdir=out_dir,
        opensn_binary=kwargs["opensn_binary"],
        processor=processor,
        n_cores=kwargs["ncores"],
        partition=kwargs["partition"],
        ngpus=kwargs["ngpus"],
        environment=kwargs["environment"],
        gpu_config=kwargs["gpu_config"]
    )
    scripts = []
    for n_nodes in nodes:
        script = create_job_script(
            engine=engine,
            keys=keys,
            n_nodes=n_nodes,
            input_script=out_dir / f"weak_scaling_{n_nodes}.py",
            prefix="weak",
            output_dir=out_dir
        )
        scripts.append(script)
    with open(out_dir / "submit_jobs.sh", "w") as launch_file:
        launch_file.write("#!/bin/bash\n\n")
        for script in scripts:
            launch_file.write(f"{commands[engine]} {script.name}\n")

    # generate the Python script
    print("Generating weak scaling Python script.")
    for n_nodes in nodes:
        script_content = py_template.render(**keys, mesh_file=f"weak_scaling_{n_nodes}.msh")
        script_name = out_dir / f"weak_scaling_{n_nodes}.py"
        with open(script_name, "w") as script_file:
            script_file.write(script_content)
