# Holdout chirality and interpolation audit

All **140 reactions / 166 minimum-event candidates** completed under a
300-second per-case watchdog. Every returned mapping passed an independent
saved-family membership check and preserved its candidate's concrete signed
broken/formed bonds. No AAM search or event decoding was rerun.

[Download and open the offline 3D viewer](viewer.html). It reuses the original
white R/P viewer, internal-coordinate interpolation, playback, frame slider,
atom labels, and magenta clash highlighting. Each changed candidate has before
and corrected choices. Unchanged valid candidates have one choice. Search the
case list for **worse** to inspect peak-clash increases.

## Paired results

Counts below refer to candidates, not independently verified chemical pathways.
Each interpolation has 101 frames. A clash means a nonbonded distance below
0.70 times the summed covalent radii; pairs bonded at either endpoint are excluded,
as in the existing interpolation code. Peak count is the maximum over all frames.

| Measure | Before | After |
|---|---:|---:|
| Candidates with no detected clashes in any frame | 14 | 84 |
| Reactions with every candidate clash-free | 11 | 68 |

- Peak clash count decreased for **129/166** candidates, increased for **10/166**,
  and was unchanged for **27/166**.
- **154/166** mappings changed; all 166 corrected/retained mappings passed the
  family and event checks. No enforced ordinary or retained high-coordinate
  frame had an orientation violation.
- Peak-clash increases occur in cases 7, 9, 18, 28, 29, 30, 40, 41, 95.
  These remain visible; chirality consistency does not guarantee a collision-free
  interpolation. Endpoint posing, torsions, and the interpolation construction
  are additional factors. No energy/path optimization was performed.

## Timing and selection policy

Chirality selection consumed **200.85 CPU seconds** in total:
mean 1.210, median 0.495,
95th percentile 3.320, maximum 25.372 seconds per candidate.
Interpolation before/after consumed 406.06 CPU seconds separately.
These stage times exclude checkpoint loading, catalogue reconstruction,
independent membership checks, and HTML packaging. Four single-threaded workers
were used. There were no watchdog terminations in this final run.

Ordinary chirality and mutability query the saved family union. Additional
high-coordinate frames form a maximal feasible basis **inside the chosen saved
family** by default. `high_coordinate_scope="union"` explicitly requests the
more expensive global basis; strict high-coordinate mode also searches the union.
Single-atom reachable sets reject impossible shuffle queries before compilation;
they never authorize a mapping. Positive answers retain all correlated actions,
edge constraints, anchors, and concrete event constraints. There is no global
RMSD ranking, bijection enumeration, or additional count cap.

The initial union-wide high-coordinate experiment exposed long negative proofs
and hit watchdogs. It was stopped before the final run. This report covers the
explicit family-local policy above. The source archives contain 300 decoded
candidates in total; the 166 observed minima match the original holdout viewer's
scope. Higher-event alternatives are outside this audit, not silently discarded
by the decoder. Endpoint mappings in this holdout set are unverified; this is a
software/geometry diagnostic, not a chemical accuracy benchmark.

## Reproduction

From the repository root, with trusted saved checkpoints at `SAVED_RUN`:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  bench/experiments/chirality_holdout/run.py \
  --source SAVED_RUN --output AUDIT --minimum-only --workers 4
PYTHONPATH=src python bench/experiments/chirality_holdout/build.py \
  --repo . --source AUDIT --output VIEWER
```

`SAVED_RUN/caseN/` contains `cuts/aam.pkl.gz`, `candidates.json`, and
`result.json` from the final end-to-end run. Saved search settings were one seed,
sweep enabled, cap 2000, competition off, and event thresholds 0.5/0.3.
The decoder's existing per-case event window is retained in the evidence.

The focused suite passed **118 tests**, including tiny exhaustive checks of
ordered/correlated actions, rejection bounds, both high-coordinate scopes,
anchors, event invariance, geometry, interpolation, and API compatibility.
Full selected mappings, frame diagnostics, stage timings, and input settings are
in [`case-results.json.gz`](case-results.json.gz). Summary statistics and source
hashes are in [`summary.json`](summary.json).

Offline browser checks passed for all 140 cases and all 320 mapping views,
including scrubbing, playback, and clash highlighting, with no JavaScript errors
or external network requests. See [`browser-check.json`](browser-check.json).
