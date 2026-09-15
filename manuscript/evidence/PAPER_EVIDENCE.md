# Evidence used by the preprint

`paper_sources.json` identifies the final source, the fresh validation report and SHA-256 hashes for each numerical snapshot. The four core result snapshots come from completed rerun campaigns; additional snapshots document paired controls, comparator rescoring, and verdict audits. Figures and tables read these snapshots; building the paper does not run an experiment.

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


Golden nonrecovery audit (2026-09-13): all 11 final ten-seed misses inspected from complete saved cuts. Four extra-pair mismatches; six exclusions by relaxed orbit bounds; one joint correspondence exclusion. Original RDF references independently verified. Current-engine uncapped causal traces reproduce exclusions in cases 986/1285 and establish reference-directed feasibility. Diagnostic results do not change 1840/1851 recovery. See golden_miss_analysis.json and reports/golden_miss_analysis_20260913.

Table 3 combines the GRAFT seed/sweep configurations, SLAP default and union/sweep configurations, and all other released comparators using reference recovery over their stated outputs. Archived first-output fields are retained only for audit; no ranking is attributed to SLAP. Figure 3 reports GRAFT's 336 classes (300 from ordinary search plus 36 added by competition), its comparator coverage, deduplication, and complete-window postprocessing cost. No searches were rerun for these presentation changes.

`direction_recovery.json` adds size-ordered recovery for the default one-seed sweep and two-seed sweep. All 1,851 prepared input hashes match the campaign manifest. Size includes explicit H; equal-size endpoints use stored order for smaller-first and reversed order for larger-first. Individual verdicts are complete at these two seed settings. Source-level seed defaults are now one; measured campaigns explicitly chose seed counts and are unchanged.

`timing_comparison.json` disaggregates the existing matched Mac cohort by smaller-first, larger-first, and bidirectional policies, with mean/median/95th-percentile CPU times and recovery counts. All per-direction CPU sums reproduce the published GRAFT totals; all saved SLAP graph/mapping/export timers reproduce its per-case and cohort totals. Uncut SLAP timing uses ordinal-zero calls only. Default-comparator call timings remain a separate table with CPU-model availability stated. No searches or reference checks were rerun.

`coordinate_minima.json` derives minimum-event count relations and pattern intersections from all140complete saved decoding results and fresh SLAP comparison records. It separates minimum-to-minimum intersections(162/164sweep;140/140native on equal-count cases) from broader window coverage(166/168;155/160). Timings read saved successful search/competition,SLAPmapping/H-refinement,and previouslyauditedfullwindowdecoder records; endpointpreparation,finalsharedSLAPrescoring and aborted/replacedattempts are outside these stagecosts.
