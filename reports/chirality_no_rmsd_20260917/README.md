# Chirality selection without global RMSD

On `dev/chirality-postprocessing`, global RMSD optimization and its covariance
search machinery have been deleted. Neither the current signed-event selector
nor the analytical-family selector ranks mappings by RMSD. The former no longer
calculates or returns a rigid-fit score. The legacy R/P result retains one
fixed-correspondence fit diagnostic computed after selection.

The analytical selector prefers a feasible source witness and otherwise returns
a constraint witness. R/P processing stops at the first feasible saved
branch/event coset. Exact family, event, orientation, and anchor constraints
remain binding. See [the API contract](../../docs/CHIRALITY.md).

## Validation

105 focused tests passed. Regression checks retain an allowed source even when
another mapping has lower RMSD, prevent geometry fitting during symbolic repair
and R/P branch selection, and retain the existing family/event/chirality tests.
The self-contained public example also passed.

Eight saved/synthetic cases passed three trials each. Every returned witness
was independently verified against its saved family and concrete signed events.
These are engineering checks, not a full-dataset benchmark. Gold uses two saved
published families; Golden uses two candidates from its complete saved baseline
catalogue. The 135-atom case uses its complete saved unswept catalogue.

One CPU per case; each trial creates a fresh selector workspace. Times exclude
imports, input loading, and independent membership verification. A parent process
provides a 300-second watchdog per case. No search or full decoding campaign was
rerun and no new count caps were introduced.

| Case | Atoms | Median CPU seconds | Median elapsed seconds |
|---|---:|---:|---:|
| independent_tetrahedra_1 | 5 | 0.0026 | 0.0026 |
| independent_tetrahedra_8 | 40 | 0.0213 | 0.0213 |
| independent_tetrahedra_24 | 120 | 0.0957 | 0.0957 |
| gold_65_a | 65 | 0.0778 | 0.0778 |
| gold_65_b | 65 | 0.0779 | 0.0779 |
| structure_verification_135 | 135 | 0.1127 | 0.1127 |
| golden_15_candidate_0 | 30 | 0.0007 | 0.0007 |
| golden_15_candidate_1 | 30 | 0.0006 | 0.0006 |

Raw trial mappings, orientation diagnostics, and timing scope are in
[`results.json`](results.json).
