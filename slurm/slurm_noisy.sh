#!/bin/bash
#SBATCH --job-name=af-edg-noisy
#SBATCH --output=logs/noisy_%j.out
#SBATCH --error=logs/noisy_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load matlab
matlab -batch "run('scripts/reproduce_noisy_constraints.m')"
