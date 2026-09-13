# GRAFT preprint

The manuscript presents the final algorithm: weighted continuous fragment growth, correlated automorphism actions, saturated-placement and fragment-priority branching, unordered final fragment pairs, and symbolic signed-event decoding.

- `manuscript.pdf`: compiled paper.
- `preprint.tex`: LaTeX entry point, retaining the supplied author list and preprint style.
- `includes/paper.tex`: methods, evaluation, and discussion.
- `includes/supplement.tex`: concise reproducibility details and completeness definitions.
- `figs/fig1_algorithm.{pdf,svg,png}`: the six-panel algorithm figure, including a checked weighted-graph example distinguishing matching from event symmetry.
- `figs/fig2_golden.{pdf,svg,png}`: Golden reference-family recovery and seed cost.
- `figs/fig3_coordinate.{pdf,svg,png}`: 140-case event-class coverage and final branch counts.
- `manuscript_bundle.zip`: compact paper, source, figure, and evidence bundle.

The manuscript does not include abandoned variants, debugging history, or superseded scoring comparisons. Older research artifacts remain in Git for provenance but are not included in the paper or its compact bundle.

## Results and scope

Every reported benchmark configuration has been rerun with the final source. Golden evaluates fragment search on all 1,851 records with 1/2/3/10 seed orderings; its reference-family recovery and paired SLAP comparison are in `evidence/seed_comparison.json` and `evidence/slap_sweep.json`. Final Golden recovery is 1,834/1,851 for one and two orderings, 1,837 for three, and 1,840 for ten (two unresolved); SLAP recovers 1,796 (three unresolved). The coordinate protocol additionally evaluates local competition and full final-family decoding, with one seed ordering only. The final two-seed coordinate pipeline was not measured. Timing compares only the 1,821 same-Mac completed reactions for one/two orderings and SLAP; mixed-host three/ten-order timing is not pooled.

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
