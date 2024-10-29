#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=00:30:00
#SBATCH --output="out.txt"
#SBATCH --error="err.txt"
#SBATCH --exclusive

srun opensn -i strong_scaling.lua
