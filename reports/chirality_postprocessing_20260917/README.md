# Chirality integration: focused validation

Historical record for `fbbf6df`, before removal of global RMSD optimization.
See [the subsequent validation](../chirality_no_rmsd_20260917/README.md) for the
current selection policy and timings.

Implemented on `dev/chirality-postprocessing`, separate from publication `main`.
The new `graft.chirality.select_chiral_witness` consumes current decoded event
candidates and the original compressed action programs. See the
[API contract](../../docs/CHIRALITY.md) and
[self-contained example](../../examples/chirality/example.py).

## Correctness

All 107 focused tests pass. Tests include exhaustive small-case oracles,
correlated/ordered actions, anchors, event boundaries, planar and changing
coordination, maximal high-coordinate bases, and timeout status. The original
RMSD regression was fixed by evaluating the selected mapping's final coordinate
residual rather than reporting a cancellation-prone covariance-norm difference.

Every real/stress timing witness was independently checked for membership in a
saved family and exact equality of concrete signed bond events. The two Golden
witnesses already satisfied orientation and required no symbolic solving. The
135-atom example changed 14 assignments in its first trial and reported 25
incompatible dependent high-coordinate frames as reconfigured; it does not
claim preservation of all coordination signs. Full per-frame diagnostics and
all selected mappings are in `results.json`.

## Small timing sample

Three trials per case, one process/CPU, fresh chirality workspace per trial.
Times measure the selector, excluding imports, input loading and the subsequent
independent membership recheck. Each case runs under a parent-enforced
300-second process watchdog. No new mapping-count caps or search limits were
introduced. No AAM search or full event-decoding campaign was rerun.

| Case | Atoms | Median CPU seconds | Median elapsed seconds |
|---|---:|---:|---:|
| independent_tetrahedra_1 | 5 | 0.0030 | 0.0032 |
| independent_tetrahedra_8 | 40 | 0.0214 | 0.0215 |
| independent_tetrahedra_24 | 120 | 0.0945 | 0.0946 |
| gold_65_a | 65 | 0.0773 | 0.0774 |
| gold_65_b | 65 | 0.0749 | 0.0751 |
| structure_verification_135 | 135 | 0.1195 | 0.1197 |
| golden_15_candidate_0 | 30 | 0.0009 | 0.0009 |
| golden_15_candidate_1 | 30 | 0.0008 | 0.0008 |

The 120-atom stress case contains 24 independent three-hydrogen pools,
representing `6**24 = 4,738,381,338,321,616,896` assignments before orientation
constraints. The first recorded trial used 25 solver checks and 24 local parity
refinements. No complete permutation/group list was expanded. This is a compact
symmetry stress test, not a claim of polynomial worst-case complexity.

Gold tests use the two published selected families, not the full 2,795-family
gold archive. The Golden example uses its complete saved baseline catalogue,
but selects only the two published event candidates. These eight cases are
engineering checks, not a replacement for a full dataset benchmark.

## Scope compared with the older solver

The coordinate determinant and ordinary mutable-center principle are retained.
The new stage returns a feasible orientation-consistent witness, **not** the
globally minimum-RMSD mapping. The old analytical-family/RMSD API remains
available. Higher-coordinate bases are now built from explicit source frames
inside an ordered/union relation, rather than symmetry-orbit units of a single
analytical coset; they need not select the same basis or witness as the older
representation. Details and strict alternatives are in the API contract.

Default mutable mode reports unavoidable fixed orientation changes instead of
rejecting genuine reaction inversion. Use `mode="all"` for a stricter structure
check. Undefined orientations and absent ligand connections remain unconstrained.
This is coordinate/index chirality, not a general CIP or E/Z stereochemistry
contract. No publication results or animations were changed.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  bench/experiments/chirality_postprocessing/run.py --output /tmp/graft-chirality
python examples/chirality/example.py
```

Install the package with `.[postprocessing]`. The saved 135-atom fixture is from
the already published molecule-verification example; the other real inputs are
read from their existing repository records. `validation.json` records source
hashes and the focused test command.
