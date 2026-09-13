#!/bin/bash
#SBATCH --job-name=aam-cut-verification
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:07:00
#SBATCH --array=0-23%24
#SBATCH --no-requeue
set -euo pipefail
cd /lustre/fsw/portfolios/nvr/users/yunhengz/Projects/coordinate_alignment/validation_20260912
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
exec .venv/bin/python hpc/verify_remaining.py "$SLURM_ARRAY_TASK_ID"
