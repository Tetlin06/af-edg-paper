#!/bin/bash
#SBATCH --job-name=af-edg-state
#SBATCH --output=logs/state_%j.out
#SBATCH --error=logs/state_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load matlab
matlab -batch "run('scripts/reproduce_state_cases.m')"
