import os


def run_gmsh(gmsh_binary, geo_filename, divisor):
    """Run Gmsh on the given .geo file to produce a mesh file."""
    gmsh_command = \
        f"{gmsh_binary} -3 -v 0 -setnumber divisor {divisor} -o strong_scaling.msh {geo_filename}"
    print(f"Generating mesh file with command: {gmsh_command}")
    os.system(gmsh_command)


def create_slurm_file(job_name, opensn_binary, nnodes, ntasks):
    slurm_filename = f'slurm_{nnodes}n.sh'
    print(f"Generating SLURM file {slurm_filename}")
    """Create a SLURM job file to run OpenSn."""
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}_{nnodes}n
#SBATCH --output={job_name}_{nnodes}n.out
#SBATCH --error={job_name}_{nnodes}n.err
#SBATCH --nodes={nnodes}
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --time=01:00:00

srun {opensn_binary} -i strong_scaling.lua
"""
    with open(slurm_filename, "w") as slurm_file:
        slurm_file.write(slurm_content)


def generate_files_for_scaling_study(opensn_binary, gmsh_binary, geo_filename, nnodes, ntasks,
                                     gmsh_divisor, job_name):
    # Generate the mesh
    run_gmsh(gmsh_binary, geo_filename, gmsh_divisor)

    for nodes in nnodes:
        # Create the SLURM job file
        create_slurm_file(job_name, opensn_binary, nodes, ntasks)
    print("")


# Name of the Gmsh geometry file
geo_filename = "cube.geo"

# Name of the base OpenSn input file for this study
input_filename = "strong_scaling.lua"

# Name of the job to include in job scripts
job_name = "strong_scaling"

# Path to the OpenSn binary
opensn_binary = "/path/to/opensn"

# Path to the Gmsh binary
gmsh_binary = "gmsh"

# Node counts
nnodes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]

# Tasks (MPI ranks) per node
ntasks = 64

# Gmsh divisor
# Generates an unstructured mesh with 268,938 tets. At 64 tasks per node, this comes to ~4200 cells
# per task for 1 node and ~8 cells per task at 512 nodes.
gmsh_divisor = 39

generate_files_for_scaling_study(opensn_binary, gmsh_binary, geo_filename, nnodes, ntasks,
                                 gmsh_divisor, job_name)
