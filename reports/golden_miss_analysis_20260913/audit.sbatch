#!/bin/bash
#SBATCH --job-name=aam-miss-audit
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:06:00
#SBATCH --array=0-21%11
#SBATCH --no-requeue
set -euo pipefail
TASK_ROOT=/lustre/fs1/portfolios/nvr/projects/nvr_lpr_agentic/users/yunhengz/Projects/coordinate_alignment/validation_20260912
cd "$TASK_ROOT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
cases=(7 19 590 871 986 1228 1285 1358 1377 1475 1553)
directions=(R_to_P P_to_R)
exec timeout 300 .venv/bin/python miss_audit_20260913/audit_retained.py "$TASK_ROOT" "$TASK_ROOT/miss_audit_20260913/results" "${cases[$((SLURM_ARRAY_TASK_ID / 2))]}" "${directions[$((SLURM_ARRAY_TASK_ID % 2))]}"
