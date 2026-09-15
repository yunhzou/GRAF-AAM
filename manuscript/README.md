# GRAFT preprint

The manuscript presents the final algorithm: weighted continuous fragment growth, correlated automorphism actions, saturated-placement and fragment-priority branching, unordered final fragment pairs, and symbolic signed-event decoding.

- `manuscript.pdf`: compiled paper.
- `preprint.tex`: LaTeX entry point, retaining the supplied author list and preprint style.
- `includes/paper.tex`: methods, evaluation, and discussion.
- `includes/supplement.tex`: concise reproducibility details and completeness definitions.
- `figs/fig1_algorithm.{pdf,svg,png}`: the branching algorithm overview: RDKit vector depictions of alcohol sites and methyl acetate saponification illustrate hydrogen-preserving local matches, three conditional assignment levels, automorphism compression, competition, and event decoding.
- `evidence/competitors.json`: audited Golden results for RXNMapper, LocalMapper, Chython, Indigo, RDT, and default SLAP configurations; Table 3 compares their returned-output reference recovery with all evaluated GRAFT seed/sweep configurations, without assigning a ranking to SLAP alternatives.
- `evidence/molecule_example.json`: checked SMILES, mappings, and heavy-atom event lists for the constructed illustration.
- `figs/fig2_golden.{pdf,svg,png}`: Golden reference-family recovery and direction-specific mean CPU cost.
- `figs/fig3_coordinate.{pdf,svg,png}`: GRAFT and prior-method minimum-event comparisons, minimum-count distributions, and measured stage costs on 140 cases.
- `figs/fig4_alternatives.{pdf,svg,png}`: three real GRAFT minimum-event alternatives from coordinate case 127; translucent red/blue highlights show signed WBO changes behind unchanged black bond strokes and pastel highlights show actual saved fragments. Includes a certified spectator-H shuffle. The discussion motivates allowed reaction-core assignments and mapping-dependent interpolation in a downstream TS pipeline; no geometric TS panel is shown.
- `evidence/coordinate_case127.json`: full saved-case witnesses, checked events, source fragment groups, and hashes for Figure 4.
- `evidence/multicandidate_example.json` and `evidence/output_multiplicity.json`: current-API example and archived comparator-output audit for the discussion.
- `manuscript_bundle.zip`: compact paper, source, figure, and evidence bundle.

The manuscript does not include abandoned variants, debugging history, or superseded scoring comparisons. Older research artifacts remain in Git for provenance but are not included in the paper or its compact bundle.

## Results and scope

The default configuration is one seed ordering per cut with the single-edge sweep. No-sweep configurations are ablations. Size-ordered directional recovery is derived from saved, hash-verified explicit-atom inputs and verdicts in `evidence/direction_recovery.json`; no search was rerun.

Every reported GRAFT benchmark configuration has been rerun with the final source. The completed default-comparator predictions were additionally rescored with the current strict evaluator; no competitor models were rerun for that check. Golden evaluates fragment search on all 1,851 records with 1/2/3/10 seed orderings; its reference-family recovery and paired SLAP comparison are in `evidence/seed_comparison.json` and `evidence/slap_sweep.json`. Final Golden recovery is 1,834/1,851 for one and two orderings, 1,837 for three, and 1,840 for ten (11 nonrecoveries, none unresolved); SLAP recovers 1,796 (three unresolved). The coordinate protocol additionally evaluates local competition and full final-family decoding, with one seed ordering only. The final two-seed coordinate pipeline was not measured. Timing compares only the 1,821 same-Mac completed reactions for one/two orderings and SLAP; mixed-host three/ten-order timing is not pooled.

The fresh coordinate configuration covers 166/168 minimum-event classes among returned SLAP-sweep witnesses, with every compared class recovered in 139/140 reactions. The SLAP hydrogen label families are not exhaustively decoded. Final branch grouping reduces 237,645 ordered records to 124,641 unordered fragment combinations. All 236,653 retained AAM families are completely decoded or certified within the fixed event windows. The complete final-decoding pass takes 11.60 wall minutes and 25.70 recorded CPU-minutes across three workers, including bounded continuation and certificate journaling. Search and competition are separate stages.

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

Runtime detail is in Table 4 (same-Mac mean/median/95th-percentile CPU by endpoint-size orientation, paired with recovery), Figure 2, and the archived default-comparator call-time table. `timing_comparison.json` retains per-case measurements and matched-cohort identities. The cluster default-mapper timings are separately scoped; no cross-host speed ranking is implied.

The140coordinate collection is a bond-event comparison with unverified mappings, not a mapping-accuracy benchmark. Tables5/6andFigure3 focus on minimum counts, patterns where minima agree, and measured stageCPU. GRAFT reaches lower/equal/higher counts in4/136/0cases versusSLAPsweep and15/125/0versusnativeSLAP. At equal minima it retains162/164sweep patterns and140/140native patterns. Its171own-minimum patterns are distinct from336classes in the larger decoded windows. No new searches or decoding were performed.
