# Installation and Run Instructions

## Installation Instructions

We have a set of easy instructions for users running Ubuntu 18.04+ (including WSL
with Ubuntu 18.04, as well as newer Ubuntu LTS)

- Easy scripts - [Easy Linux instructions](./install_linux_easy.md)
- Easy scripts - [Easy MacOS instructions](./install_macos_easy.md)

If the automated installation of dependencies fails, or if you don't have a system
supporting it, then follow the instructions below: 

For Linux machines - [Linux installation instructions](./install_linux.md)  
For MacOS machines - [MacOS installation instructions](./install_macos.md)

The executable is typically found in `opensn/build/lua`

## Folder Layout
The subdirectory structure is as follows:

- [doc](./doc) contains 
  - scripts to generate the code documentation (these scripts were actually used during the installation instructions, above),
  - the generated documentation;
- [external](./external) contains external library source code, such as VTK,
- [framework](./framework) contains the source code that is agnostic of the physics, such as math, mesh, MPI, ...
- [modules](./modules) contains the source code for specialized physics modules,
- [resources](./resources) contains python scripts, notes, data, and other misc. items,
- [tests](./tests) contains tests for the regression suite (these can also serve as examples/tutorials),
- [tutorials](./tutorials) contains tutorials. 


## Run Instructions

###  Serial run (interactive)
Simply type:
```bash 
path_to_exec/opensn -i path_to_input_file/input_filename.lua
```
where 
- ```path_to_exec``` is the path to the directory where the executable was created and saved,
- ```opensn``` is the (default) filename for the executable,
- ```path_to_input_file``` is the path to the directory where the input file is located,
- ```input_filename.lua``` is the name of the LUA input file

It is also possible to run in batch mode on a cluster or a supercomputer. See below.

### Parallel run (interactive)
In interactive mode (i.e., on your own machine), simply type
```bash 
mpiexec -n N path_to_exec/opensn -i path_to_input_file/input_filename.lua
```
where ```N``` is the requested number of parallel processes.

### Parallel run (batch)

In batch mode, you have to use the job submission of the cluster you intend 
on using.

For example, a slurm submission script may look like this:
```bash
#!/usr/bin/bash
#
#SBATCH -J {file_name} # Job name
#SBATCH -o tests/{file_name}.o # output file
#SBATCH -e tests/{file_name}.e # error file
#SBATCH -p skx-normal_ # Queue (partition) name
#SBATCH -N {num_procs // 48 + 1} # Total # of nodes
#SBATCH -n {num_procs} # Total # of mpi tasks
#SBATCH -t 00:05:00 # Runtime (hh:mm:ss)
#SBATCH -A Massively-Parallel-R # Allocation name (req'd if you have more than 1)
```
