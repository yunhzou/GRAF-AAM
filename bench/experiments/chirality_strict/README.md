# Strict chirality audit against the archived 140-case viewer

This experiment uses the existing final AAM archives and decoded minimum-event candidates. It does not rerun search or decoding. It compares the uploaded old mappings against an independently implemented orientation test, certifies new witnesses in the saved AAM families, and queries each same-event old mapping that passes the strict geometry test.

The report is in `reports/chirality_strict_audit_20260917/`. The old viewer is a reference, not ground truth. Counts of defined orientation reversals are not counts of chemically invalid reactions.

## Snapshot inputs

Set `GRAFT_AUDIT_WORKSPACE` to a workspace containing these existing artifacts:

- `outputs/final-end-to-end-optimized-20260916/caseN/{cuts/aam.pkl.gz,candidates.json,result.json}`: trusted saved archives and decoded results.
- `work/chirality-comparison/old-<reaction name>.json`: per-reaction `DATA` objects extracted from the old uploaded viewer (its `CASES` array contains gzip/base64 HTML payloads). Parse these as data; do not execute scripts.
- For the comparison viewer only, `work/full140_inputs/N/input.json` and `outputs/chirality-holdout-minimum-20260917/caseN/viewer-data.json.gz`: endpoints and previously computed raw interpolations.

These large private inputs are not added to Git. The report's compressed case records contain selected mappings, orientation diagnostics, and exact membership results. The offline HTML embeds the endpoint coordinates and displayed trajectories.

```sh
export GRAFT_AUDIT_WORKSPACE=/absolute/path/to/workspace
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PYTHONPATH=src python bench/experiments/chirality_strict/run.py
PYTHONPATH=src python bench/experiments/chirality_strict/build_viewer.py
```

The default run uses four isolated workers, 300-second per-case watchdogs, and 240-second soft selector watchdogs. There are no new branch or solution caps. Inconclusive cases must remain `unknown`/unfinished and must not be counted as infeasible. To audit selected cases use `--cases 3 31 68 95`. Use a fresh output directory (or relocate the old audit folder) after changing code or inputs: completed case files are reused.

The independent geometry oracle directly computes determinants and normalized signs, without the production solver's orientation helpers. The positive witness is separately checked with `query_path` fixing every atom and recomputing concrete signed events. Small exhaustive tests in `tests/test_chirality_postprocessing.py` validate negative answers and the rejection bounds; four real-case full-solver pilots agree with the accelerated negative answers.

Run the existing offline browser test with `bench/experiments/chirality_holdout/check_viewer.cjs`. Raw mappings displayed for infeasible candidates are explicitly labeled **NO STRICT SOLUTION · raw display only**. They are not successful corrections.
