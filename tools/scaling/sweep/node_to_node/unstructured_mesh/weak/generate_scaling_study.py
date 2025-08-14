import os


def run_gmsh(gmsh_binary, geo_filename, divisor, nnodes, ntasks):
    """Run Gmsh on the given .geo file to produce a mesh file."""
    gmsh_command = f"{gmsh_binary} -3 -v 0 -setnumber divisor {divisor} " \
        f"-o {nnodes}n_{ntasks}rpn.msh {geo_filename}"
    print(f"Generating mesh file with command: {gmsh_command}")
    os.system(gmsh_command)


def update_input_file(input_filename, nnodes, ntasks):
    """Update the filename in the input file and create a new input file for each node count."""
    updated_filename = f"weak_scaling_{nnodes}n.py"
    print(f"Generating input file {updated_filename}")

    with open(input_filename, "r") as file:
        content = file.readlines()

    with open(updated_filename, "w") as file:
        for line in content:
            if line.strip().startswith('FromFileMeshGenerator(filename="meshfile")'):
                file.write(f'        FromFileMeshGenerator(filename="{nnodes}n_{ntasks}rpn.msh")\n')
            else:
                file.write(line)


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

srun {opensn_binary} -i weak_scaling_{nnodes}n.py
"""
    with open(slurm_filename, "w") as slurm_file:
        slurm_file.write(slurm_content)


def generate_files_for_scaling_study(opensn_binary, gmsh_binary, input_filename, geo_filename,
                                     nnodes, ntasks, gmsh_divisors, job_name):
    for nodes, divisor in zip(nnodes, gmsh_divisors):
        print(f"Generating inputs for {nodes} nodes...")

        # Generate the mesh
        run_gmsh(gmsh_binary, geo_filename, divisor, nodes, ntasks)

        # Create input file for each run
        update_input_file(input_filename, nodes, ntasks)

        # Create the SLURM job file
        create_slurm_file(job_name, opensn_binary, nodes, ntasks)

        print("")


# Name of the Gmsh geometry file
geo_filename = "cube.geo"

# Name of the base OpenSn input file for this study
input_filename = "weak_scaling.py"

# Name of the job to include in job scripts
job_name = "weak_scaling"

# Path to the OpenSn binary
opensn_binary = "/path/to/opensn"

# Path to the Gmsh binary
gmsh_binary = "/path/to/gmsh"

# Node counts
nnodes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]

# Tasks (MPI ranks) per node
ntasks = 64

# Divisors for creating mesh files. These values go together with node counts
# and tasks per node to generate meshes with the appropriate number of cells
# per task.
# 64 tasks per node with 256 cells per task
gmsh_divisors = [15, 19, 24, 31, 39, 49, 62, 78, 98, 124]
# 96 tasks per node with 256 cell per task
# gmsh_divisors = [17, 22, 28, 35, 44, 56, 71, 89, 113, 137]
# 64 tasks per node with 2048 cells per task
# gmsh_divisors = [30, 39, 48, 62, 78, 98, 124]

generate_files_for_scaling_study(opensn_binary, gmsh_binary, input_filename, geo_filename,
                                 nnodes, ntasks, gmsh_divisors, job_name)
