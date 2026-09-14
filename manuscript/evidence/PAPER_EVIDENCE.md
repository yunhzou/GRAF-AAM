# Evidence used by the preprint

`paper_sources.json` identifies the final source, the fresh validation report and SHA-256 hashes for each numerical snapshot. All four snapshots come from completed rerun campaigns. Figures and tables read these snapshots; building the paper does not run an experiment.

| Snapshot | Manuscript use |
|---|---|
| `seed_comparison.json` | All 1,851 Golden records at each seed setting; verified outcomes; same-Mac workflow CPU for one/two orderings and SLAP on 1,821 completed reactions |
| `slap_sweep.json` | Fresh pinned Golden SLAP sweep, including paired recovery differences |
| `competition_final.json` | Fresh 140-case search/competition, complete AAM event windows and fresh SLAP witness comparison |
| `final_dedup.json` | Fresh ordered/unordered counts and complete final decoding of all 140 cases, with attempt-inclusive timing |

Figure 1 uses RDKit vector molecule depictions. The builder checks the illustrated partial mappings, matching-response labels, and signed heavy-atom events. Figures 2 and 3 use the measured snapshots above.

Golden measures the fragment-search stage. The coordinate protocol additionally includes competition and final decoding. The coordinate collection has no annotated mappings and was used during development, so comparator event-class coverage is not chemical accuracy. Its SLAP inputs restore and verify the original XYZ-derived graph adjacency; XYZ preparation itself is not rerun. The comparator classes come from returned H-refined witnesses, not exhaustive decoding of its hydrogen label families.

A verified witness establishes recovery. Unresolved computation is not a negative result. Complete decoding certifies only the saved AAM families within the reported event window. Timings describe their recorded stage and resource scope and do not imply an isolated end-to-end speed comparison.

The main text excludes development history. The full validation report retains execution limits, continuation attempts, source hashes and the explanation of tied SLAP hydrogen-refinement witnesses. Earlier research records remain outside this paper bundle.

The main algorithm figure uses a constructed methyl acetate hydrolysis example, independently checked in `molecule_example.json`. It is not a benchmark trajectory. Formal weights and heavy-atom event counts are illustrative; the reported benchmark evaluates explicit hydrogens.

`competitors.json` restores the completed seven-configuration Golden reproduction. All 12,957 saved records were freshly rescored with the current evaluator, with unchanged outcomes. These default explicit-output results have a different scope from the final GRAFT family search and SLAP sweep. Full audit and per-case checks: `reports/golden_competitor_recheck_20260913/`.

Paired no-sweep controls are included in the main Golden results table. Evidence: `reports/golden_unswept_20260913/`; GRAFT 1,489/1,851, SLAP 1,661/1,851 (one unresolved). No new mapping searches were run. The molecular figure now checks alcohol-site matches including attached H and shows three conditional assignment levels.


2026-09-13 complete-cut audit: ten-seed Golden recovery is 1,840/1,851, with 11 reference mappings absent and zero unresolved. Cases 590 and 1358 have complete negative cut coverage. Original failed monolithic execution records remain in the campaign report. See checkpoint_verdict_audit.json and reports/checkpoint_memory_fix_20260913; no search reruns or change in recovery percentage.
