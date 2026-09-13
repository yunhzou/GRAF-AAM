#!/bin/bash
#SBATCH --job-name=aam-summary
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
cd /lustre/fsw/portfolios/nvr/users/yunhengz/Projects/coordinate_alignment/validation_20260912
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec .venv/bin/python hpc/aggregate.py
