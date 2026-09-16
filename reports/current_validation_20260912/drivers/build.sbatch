#!/bin/bash
#SBATCH --job-name=aam-build-check
#SBATCH --partition=cpu_short
#SBATCH --account=nvr_lpr_agentic
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
cd /lustre/fsw/portfolios/nvr/users/yunhengz/Projects/coordinate_alignment/validation_20260912
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
AAM_PY="$PWD/.venv/bin/python"
cd engine
"$AAM_PY" setup.py build_ext --inplace
"$AAM_PY" native/build_engine.py
cd ..
cp engine/src/rxn_core/*.so optimized/src/rxn_core/
PYTHONPATH=optimized/src "$AAM_PY" -m pytest optimized/tests/test_action_reduction.py optimized/tests/test_final_branches.py optimized/tests/test_family_query.py optimized/tests/test_event_certificates.py optimized/tests/test_event_patterns.py -q
PYTHONPATH=engine/src "$AAM_PY" -m pytest engine/tests/test_native.py engine/tests/test_native_group_ops.py -q
"$AAM_PY" score_checkpoints.py test
"$AAM_PY" hpc/build_proof.py
