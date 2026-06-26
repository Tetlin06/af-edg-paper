#!/bin/bash
#SBATCH --job-name=af-edg-exact
#SBATCH --output=logs/exact_%j.out
#SBATCH --error=logs/exact_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load matlab
matlab -batch "run('scripts/reproduce_exact_constraints.m')"
