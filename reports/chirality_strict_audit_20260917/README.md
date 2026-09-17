# Strict chirality audit: 140 holdout reactions

The implementation now reports strict orientation feasibility over the complete saved AAM family union. It does not silently relax a contradictory frame or restrict the default search to whichever family supplied the first witness. Exact and partial membership queries use the same strict predicate. This work remains on the experimental chirality branch; main was not changed.

[Open the offline 3D comparison viewer](viewer.html). Old mappings appear as reference displays with their independently counted orientation reversals. Current mappings labeled **STRICT VERIFIED** passed both independent geometric validation and exact family/event membership checks. Entries marked **NO STRICT SOLUTION · raw display only** show an uncorrected mapping solely for inspection.

## Results

The 166 current minimum-event candidates represent 140 reactions. For each candidate, the search fixes its concrete broken and formed edges and considers all current saved complete AAM families.

| Strict result | Candidates |
|---|---:|
| At least one fully preserving saved mapping | 96 |
| Proven no fully preserving saved mapping with those concrete events | 70 |
| Unresolved or timed out | 0 |

At the reaction level, all tested minimum-event candidates pass in 81 reactions, some pass in 7, and none pass in 52. The last number does **not** establish that no higher-event candidate, unsearched mapping, or physically valid reaction pathway exists. AAM search itself was not repeated or made exhaustive; negatives concern its saved families with the specified concrete event pattern.

The 96 accepted mappings were independently checked for complete bijection, element compatibility, exact family membership, unchanged concrete signed events, and every defined requested orientation. All passed. All 70 negative results visited or soundly ruled out every saved family; no timeout was counted as a negative.

## Checking the old viewer

Under the same full geometric rule, 93 old mappings preserve every defined orientation; 73 reverse at least one. The old metadata's zero-violation count excluded deliberately reconfigured higher-coordinate frames and fixed ordinary inversions. Thus the old viewer was not a certificate of full preservation.

Among the 93 geometrically preserving old mappings, 82 have a matching current concrete-event candidate. Exact post-chirality queries accept 76; 6 are already absent from the current saved AAM family union. The remaining 11 have no matching current concrete-event candidate and are outside that paired query. Every old mapping that passes geometry and current AAM membership is retained by the corrected strict filter. Choosing another witness does not remove it.

The six strictly preserving old mappings absent from the current saved relation are cases 11, 29, 30, 96, 111, and 119. This is an archive/search-family difference, not evidence that chirality rejected those mappings. We did not enlarge the AAM relation or import an old mapping as an unauthorized solution. The uploaded HTML does not contain the full old search archive/settings needed to attribute this difference conclusively.

Earlier, eight old mappings were found to violate the new run's chosen high-coordinate frame basis. That alone did not prove the old choices correct: full checks show that they involved orientation reconfiguration. The actual software problems were a witness-dependent default relaxation, restricting subsequent selection to one family, presenting relaxed feasibility as successful preservation, and activating some extra frames from an unrelated atom of the same element. Those behaviors have been corrected or made explicitly opt-in and labeled `relaxed`.

## Interpretation

This is a **coordinate/index-orientation contract**, including explicit chemically equivalent atoms. It is not a CIP or E/Z assignment. Three persistent neighbors define a center-relative signed volume; four define an affine ligand tetrahedron. Centers high-coordinate at either actual mapped endpoint additionally require applicable persistent three/four-neighbor signs to agree. Undefined endpoint frames and lost connections are inactive.

Real reactions may invert a center or rearrange a metal coordination environment. Therefore “no strict solution” is not a verdict that a reaction or old mapping is chemically wrong. It means the specified full-preservation constraints cannot all be met within the saved relation. Explicit historical relaxation remains available, reports every full-orientation violation, and returns `relaxed` rather than falsely claiming a strict correction.

Strict chirality is also not a clash optimizer: 53 of the 96 certified witnesses have no clashes anywhere in the 101-frame interpolation; the other 43 still have overlaps. Those overlaps must not be relabeled as chirality violations or hidden by dropping difficult candidates. The old and new complete candidate collections differ, so these counts are not an accuracy ranking between methods.

## Runtime and verification

Strict selection took **509.83 CPU seconds total**, **3.07 CPU seconds per candidate on average**, with a **0.258-second median**. This includes proving the negative cases and excludes loading checkpoints, auditing old mappings, independent validation, interpolation, and viewer generation. It is not directly comparable to the earlier permissive selector's timing.

Four workers and 300-second per-case watchdogs were used, with a 240-second soft budget per selector. There were no unresolved cases. No branch or solution cap was added; no AAM search or decoding was rerun. Rejection bounds only reject families whose singleton reachable images force a reversal. Positive answers always use the full correlated solver. Exhausted strict-search models are released.

Four real-case pilots without the forced-orientation rejection bound (3, 31, 68, 95) agree with the accelerated infeasibility results. **103 focused tests passed**, and the public example plus a relocated benchmark-runner smoke test executed successfully. The focused suite covers small exhaustive action-program oracles, preferred and external witnesses, cross-family subset membership, partial assignments, anchors, exact events, ordinary/fixed inversions, higher coordination, degenerate geometry, reference-dependent relaxation labels, timeouts, and bound equivalence.

## Reproduction and data

See [the benchmark runner](../../bench/experiments/chirality_strict/README.md) for snapshot inputs and commands. `case-results.json.gz` stores all 140 audited case records, old geometric diagnostics, exact old-witness queries, selected mappings, feasibility statuses, family certificates, and timing. `summary.json` stores counts and source hashes; `viewer-stats.json` stores interpolation clashes, separated from feasibility. The report uses no global RMSD ranking and changes no endpoint coordinates.

## Cases without a strict solution among tested minima

3, 5, 6, 7, 8, 9, 15, 17, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 44, 45, 46, 50, 51, 52, 53, 54, 55, 56, 57, 68, 69, 80, 86, 90, 93, 94, 95, 97, 99, 100, 101, 109, 110, 111, 114, 115, 116, 117, 123, 129, 132, 136.

These are zero-based dataset indices, also shown in the viewer.

Offline browser validation exercised all 140 cases and 332 mapping views, including playback, scrubbing, and clash highlighting: no JavaScript errors and no external requests. Catalog labels distinguish mapping views from unique event candidates and do not display the old misleading zero-violation total.
