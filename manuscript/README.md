# GRAFT preprint

The manuscript presents the final algorithm: weighted continuous fragment growth, correlated automorphism actions, saturated-placement and fragment-priority branching, unordered final fragment pairs, and symbolic signed-event decoding.

- `manuscript.pdf`: compiled paper.
- `preprint.tex`: LaTeX entry point, retaining the supplied author list and preprint style.
- `includes/paper.tex`: methods, evaluation, and discussion.
- `includes/supplement.tex`: evaluation scope, completeness argument, and supplementary results.
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md): full technical settings and links to code, input/source pins, and per-case evidence.
- `figs/fig1_algorithm.{pdf,svg,png}`: the branching algorithm overview: RDKit vector depictions of alcohol sites and methyl acetate saponification illustrate hydrogen-preserving local matches, three conditional assignment levels, automorphism compression, cut sweep, and event decoding.
- `evidence/competitors.json`: audited Golden results for RXNMapper, LocalMapper, Chython, Indigo, RDT, and default SLAP configurations; Table 3 compares their returned-output reference recovery with all evaluated GRAFT seed/sweep configurations, without assigning a ranking to SLAP alternatives.
- `evidence/molecule_example.json`: checked SMILES, mappings, and heavy-atom event lists for the constructed illustration.
- `figs/fig2_golden.{pdf,svg,png}`: Golden reference-family recovery and direction-specific mean CPU cost.
- `figs/fig3_coordinate.{pdf,svg,png}`: GRAFT and prior-method minimum-event comparisons, minimum-count distributions, and measured stage costs on 140 cases.
- `figs/fig4_alternatives.{pdf,svg,png}`: Golden case 9, with the recovered reference and two alternatives. Red × marks and green inward arrows identify bond changes; pastel colors identify actual saved fragments. The discussion motivates allowed reaction-core assignments and mapping-dependent endpoint interpolation as downstream uses.
- `evidence/golden_case9.json` and `golden_case9_verification.json`: full witnesses, source fragment groups, source hashes, original-RDF reference verification, and H-event lower bounds for the selected heavy mappings.
- `evidence/multicandidate_example.json` and `evidence/output_multiplicity.json`: current-API example and archived comparator-output audit for the discussion.
- `manuscript_bundle.zip`: compact paper, source, figure, and evidence bundle.

The manuscript does not include abandoned variants, debugging history, or superseded scoring comparisons. Older research artifacts remain in Git for provenance but are not included in the paper or its compact bundle.

## Results and scope

`evidence/paper_sources.json` records the source commit, snapshot hashes and full report location. The figure and table builder requires complete fresh-campaign evidence; pending outcomes cannot silently enter a final table.

## Build

Install Python 3.12, Tectonic, and the packages in `requirements-build.txt` in an isolated environment. Then:

```sh
MANUSCRIPT_PYTHON=/path/to/python bash scripts/build.sh
```

The build regenerates figures and tables from checked-in evidence, compiles the paper, verifies numerical claims and references, and creates the compact bundle. It performs no mapping searches. The first Tectonic build may download TeX packages. A saved `preprint.bbl` accompanies the source bundle.

For visual review, render the PDF with Poppler:

```sh
pdftoppm -scale-to 1400 -png manuscript.pdf build/page
```

`EDITORIAL_NOTES.md` records the remaining author-supplied metadata and dataset provenance needed for submission. These are not fabricated in the paper. This is an arXiv-oriented manuscript, not a submitted arXiv record.

Paired no-sweep controls are included in the main Golden results table. Evidence: `reports/golden_unswept_20260913/`; GRAFT 1,489/1,851, SLAP 1,661/1,851 (one unresolved). No new mapping searches were run. The molecular figure now checks alcohol-site matches including attached H and shows three conditional assignment levels.

Runtime detail is in Table 4 (same-Mac mean/median/95th-percentile CPU by endpoint-size orientation, paired with recovery), Figure 2, and the archived default-comparator call-time table. `timing_comparison.json` retains per-case measurements and matched-cohort identities. Default-mapper timings are reported separately because their CPU model is not recorded in the archived summaries; no cross-CPU speed ranking is implied.

The 140-reaction coordinate collection compares bond events under unverified mappings; it is not a mapping-accuracy benchmark. The cap control returns complete mappings for 138 reactions at cap 100 and all 140 at cap 2,000; cases 123 and 125 are the two search-completion losses. Table 5 reports the baseline-only bond-event comparison at cap 2,000, Table 6 reports archived stage costs, and Table A5 reports fresh full-pipeline timing. Against SLAP sweep and native SLAP, lower/equal/higher event counts are 4/136/0 and 15/125/0, respectively. Final event-pattern evidence is in `reports/published_baseline_20260915/`; older competition-enabled decoding reports are excluded.

Figure 4 uses Golden case 9: the recovered reference and two alternative mappings. Compact red × marks for losses and paired green inward arrows for gains sit on the affected bonds in Figures 1 and 4, retaining fragment colors and the original black bond lines. Golden labels concern heavy atoms; totals explicitly include unannotated H assignments. The three selected witnesses are not an exhaustive decoded set. Evidence: `golden_case9.json` and `golden_case9_verification.json`.

## 3D research preview

[Grow. Branch. Decode.](animations/graft_research_preview/index.html) is a 22-second baseline-only 3D film of Golden case 15. Conditional fragment placements grow into a tree and the two shown branches decode to two distinct event classes. The complete baseline catalogue has 44 unordered branches, 96 families and nine classes through six events; two classes are shown. [MP4](animations/graft_research_preview/graft-grow-branch-decode.mp4) · [GIF](animations/graft_research_preview/graft-grow-branch-decode.gif).

The completed same-CPU Golden cap ablation is in the appendix (`evidence/golden_cap_ablation.json`). It includes cap 100 and cap 2,000 accuracy, unknown outcomes, and search CPU for all four swept seed settings, both uncut controls, and SLAP baselines. All paired timings use the same 1,403 reactions; all-attempt recorded cost is shown separately. Full results are in `reports/golden_controlled_20260915/final/`. The 3D film uses the verified baseline Golden example in `reports/published_baseline_20260915/film/`.

## Published pipeline

The published method uses one-seed cut sweep followed by separate decoding. Competition is excluded. Golden results are unchanged. Baseline-only coordinate evidence and decoding checks are in [the publication report](../reports/published_baseline_20260915/). At cap 2,000, GRAFT recovers 162/168 SLAP-sweep patterns and 151/160 native-SLAP patterns within the fixed windows. At equal minima, the corresponding counts are 158/164 and 138/140. The output has 300 event patterns, including 166 at GRAFT's own minima.

Final end-to-end timing is recorded in `evidence/end_to_end_timing.json` and Table A5: all 140 reactions complete within the five-minute watchdog, preserving all 300 archived event classes and all 120,052 saved families. Mean/median CPU seconds are 1.257/0.300 for search, 4.647/0.288 for decoding, and 6.073/0.812 end to end. These are per-reaction CPU times, not divided by the four workers. Golden timings retain their explicitly labeled search-only scope. The final paper results, evidence, and bundle are ready; no further run is pending for these reported claims.
