# rxn_core

The default GRAFT pipeline uses one seed ordering per cut, the uncut plus single-edge sweep, and branch cap 100. Match with `search_aam`, then decode bond-event alternatives separately. See the [Python API example](docs/PYTHON_API.md). Fragment competition is an experimental, explicit opt-in extension and is off by default; it is not part of the published method.

Symmetry-aware WBO atom mapping, analytical R/P alignment, and
mechanism-local transition-state analysis.

Paper benchmarks: [results and artifact index](reports/README.md) ·
[benchmark scripts and reproduction guide](bench/README.md).

Interactive views use [two shared styles](docs/VIEWERS.md): the original white
R/P/TS presentation and the catalog/results presentation.

## Grow, branch, decode

![GRAFT: Golden fragment growth, branching and two distinct decoded event classes](manuscript/animations/graft_research_preview/graft-grow-branch-decode.gif)

A 22-second Golden example: fragments grow into conditional branches, and the two grown branches decode to two distinct event classes on the right. Two decoded event classes are illustrated, including the reference-equivalent class; the saved catalogue contains nine. Red × marks show breaking or weakening; green inward arrows show forming or strengthening.

[Full-resolution MP4](manuscript/animations/graft_research_preview/graft-grow-branch-decode.mp4) · [Offline interactive viewer](manuscript/animations/graft_research_preview/index.html) · [Example and provenance](manuscript/animations/graft_research_preview/README.md)

## Current work and stable versions

See the [branch guide](docs/BRANCHES.md) for preserved pre-acceleration baselines,
active development, and archived experiments. The
[paper branch](https://github.com/yunhzou/coordinate_alignment/tree/paper/continuous-fragment-growth/manuscript)
contains the manuscript PDF, figures, animations, and reproducible figure data.

## Design

The computational API consists of immutable typed stages:

```text
search_aam
    -> group_mechanisms                 # optional
        -> compile_mechanism_families
        -> select_rp_mappings
            -> analyze_transition_state
```

AAM is the authoritative source of mapping information. Its result retains
a fragment-decision graph with shared prefixes/reconvergence, compressed
placements, exact target generators, seed/cut provenance, and cap records.
Mechanisms are an optional post-processing result, not the core container.
`match_fragment` is independently reusable by AAM and retro detection.
R/P and TS processing consume those objects;
they do not reconstruct an alternative AAM model from serialized records.

See [the search-graph API](docs/AAM_SEARCH_GRAPH_API.md) for the object model,
conditional fragment API, persistence, and path replay;
[ALGORITHM.md](ALGORITHM.md) covers downstream algorithms.
See [retro detection and assembly](docs/RETRO_ASSEMBLY.md) for geometric
building-block recommendation using the same matcher and saved AAM graphs.
The separate [big-block / gap-first beta](docs/RETRO_BETA.md) defers augmentation
until a reactant is selected; it does not replace the full workflow.

## Install

```bash
python -m pip install -e .
```

Installation builds the native bookkeeping extensions and requires a C++17
compiler; build isolation installs pybind11. Python still owns the AAM search
graph, fragment hierarchy and retro pipeline. The `_group_ops` extension
accelerates exact permutation/occupation operations, not search-policy changes.

Install xTB separately only when endpoint WBO matrices must be computed. The
typed core API accepts already materialized coordinates and WBO matrices and
does not invoke xTB.

## Python API

Start with the short [AAM notebook](docs/AAM_SIMPLE.ipynb): matching, unique bond-event candidates, certified symmetry queries, and py3Dmol inspection. The [Python API guide](docs/PYTHON_API.md) lists anchors, directions, conditional matching, all configuration controls, and current chirality limitations. Install its dependencies with `python -m pip install -e ".[notebook]"`.

```python
from rxn_core import (
    AAMProblem,
    AAMSearchConfig,
    MolecularEndpoint,
    TransitionStateTarget,
    VibrationalModes,
    analyze_transition_state,
    group_mechanisms,
    compile_mechanism_families,
    search_aam,
    select_rp_mappings,
)

reactant = MolecularEndpoint(elements_R, xyz_R, wbo_R, label="R")
product = MolecularEndpoint(elements_P, xyz_P, wbo_P, label="P")
problem = AAMProblem(reactant, product, name="reaction")

aam = search_aam(problem, AAMSearchConfig(), workers=8,
                 intermediate_dir="alignment/aam_search")
grouped = group_mechanisms(aam)
families = compile_mechanism_families(
    grouped, workers=8, minimum_events_only=True)
rp = select_rp_mappings(families)

target = TransitionStateTarget(
    MolecularEndpoint(elements_TS, xyz_TS, wbo_TS, label="TS"),
    VibrationalModes(frequencies, normal_modes),
)
ts = analyze_transition_state(rp, target)
```

For R/P only, `align_reaction(problem, workers=8)` is the convenience
composition of the search, grouping, family, and selection stages. They remain available
when callers need to inspect, cache, audit, or transform AAM information.

Serialization, CLI workflows, and self-contained HTML views are typed
artifact adapters, not computational data models.

## CLI

NPZ endpoint files contain `elements`, `coordinates`, and `wbo` arrays:

```bash
rxn-core \
  --stage rp \
  --reactant-npz R.npz \
  --product-npz P.npz \
  --workers 48 \
  --post-workers 8 \
  --output alignment
```

Existing xTB cache directories containing one XYZ and a `wbo` file can be
used directly:

```bash
rxn-core --reactant-cache cache/R --product-cache cache/P \
  --workers 48 --output alignment
```

The output contains `rp.json`, one `R.xyz` and `P_aligned.xyz` pair per
mechanism, a reusable `reaction.json`, and a self-contained `view.html`.
The TS verifier/scorer can then be entered independently without R/P search:

```bash
rxn-core --stage ts \
  --reactant-npz R.npz --product-npz P.npz \
  --reaction-json alignment/reaction.json \
  --target-npz guess_1.npz --target-npz guess_2.npz \
  --output ts_scores
```

Use `--stage full` to compose both stages in one process. Target NPZ files
contain `elements`, `coordinates`, `wbo`, `frequencies`, and `modes`.

The executed [tutorial notebook](docs/TUTORIAL.ipynb) demonstrates the staged
Python API, CPU controls, AAM inspection, R/P alignment, and TS scoring.

## Result hierarchy

```text
AAMResult
`- AAMSearchGraph
   |- contexts, states, transitions, and stop/cap evidence
   `- compressed placements, conditioned generators, and provenance

DecodedEvents                              # separate signed-event post-processing
`- EventCandidate[]                        # one witness per unique event class
   |- event edges, count, and supporting family IDs
   `- symmetry evidence and conditional shuffle queries

MechanismResult / AnalyticalAAMResult       # optional geometry/TS pipeline
`- analytical mapping families per mechanism

RPResult
`- selected mapping per mechanism, with index-chirality diagnostics

TSResult
`- R->TS and P->TS CoreAAMResult plus scored core tuples
```

## Tests

```bash
.venv/bin/pytest -q
```

The suite includes a non-empty TS integration case that performs endpoint
AAM, analytical-family compilation, R/P selection, two partial core searches,
endpoint-consensus merging, and imaginary-mode scoring.
