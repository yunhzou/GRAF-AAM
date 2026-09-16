#!/bin/bash
#SBATCH --job-name=aam-miss-trace
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --array=0-1%2
#SBATCH --no-requeue
set -euo pipefail
TASK_ROOT=/lustre/fs1/portfolios/nvr/projects/nvr_lpr_agentic/users/yunhengz/Projects/coordinate_alignment/validation_20260912
cd "$TASK_ROOT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
cases=(986 1285)
exec timeout 240 .venv/bin/python miss_audit_20260913/trace_causes.py "$TASK_ROOT" "$TASK_ROOT/miss_audit_20260913/traces" "${cases[$SLURM_ARRAY_TASK_ID]}"
