# GRAFT — Grow. Branch. Decode.

A 48-second, silent 3D research preview. Open `index.html` for the offline,
scrubbable film, or play `graft-grow-branch-decode.mp4`. The looping GIF is a
smaller preview; the MP4 preserves label clarity and smooth camera motion.

The film uses the repository's bundled 3Dmol renderer, black endpoint bonds,
matched fragment colors, and a moving cue connecting an actual atom pair.
The tree progressively reveals both an initial fork and a conditional second
fork. The closing cards show complete matched endpoint pairs and separate
bond-event labels, without recoloring the fragment bonds.

## Storyboard

| Time | Scene |
| --- | --- |
| 0–3 s | Introduce the endpoint pair |
| 3–7 s | Grow the first fragment under conditional isomorphism |
| 7–10 s | Retain two placements of that fragment |
| 10–19 s | Grow the next fragment under the first placement |
| 19–23 s | Reveal a second, conditional fork; preview each branch |
| 23–28.5 s | Complete the retained fragment paths |
| 28.5–33 s | Extract recorded witnesses from the retained representation |
| 33–48 s | Reveal three distinct bond-event patterns and their mappings |

## Scientific scope

The source is `reports/pr7_search_trajectory_20260911/trace.json`, coordinate
case 101 (PR7). This is an archived visualization example, **not the ongoing
cap-2,000 Golden experiment**. It uses one saved search context from seed order
15 of a 30-order run, branch cap 1,000, matching tolerance 1.0, event tolerance
0.5, graph floor 0.2, and the C18–O21 sweep constraint. Atom indices are the
original zero-based source identities; the product is labeled with the
corresponding source identities.

The 12 nodes and 11 edges are copied from the saved graph. Its three terminal
witnesses each have four events, but their atom-indexed event patterns differ.
They are selected recorded witnesses, **not an exhaustive symmetry decode**.
The animation extracts each saved terminal correspondence; it does not run a
new search or benchmark the decoder. Temporary growth candidates can contain
symmetry blocks. Those blocks are retained in `film-data.json` and are not
turned into fictitious structural tree branches.

Colors encode each witness's matched fragment sets. Product coordinates are
rigidly oriented for visibility; the slow camera movement is purely visual.
No endpoint interpolation, molecular dynamics, transition state, or confirmed
reaction pathway is claimed. The cut is a search condition, not a declaration
that the corresponding bond must break.

`science-validation.json` records independent WBO-event recomputation and
mapping checks. `browser-validation.json` records timeline/browser checks.
`film-data.json` contains full mappings, events, retained graph, growth frames,
source hash, engine commit, and archive provenance.

## Rebuild

From the repository root, with NumPy, Playwright, Chrome/Chromium and FFmpeg
installed:

```sh
python manuscript/scripts/research_film/build_film.py \
  --repo . --output manuscript/animations/graft_research_preview
node manuscript/scripts/research_film/render_film.cjs \
  manuscript/animations/graft_research_preview video
python manuscript/scripts/research_film/export_gif.py \
  manuscript/animations/graft_research_preview
```

`PLAYWRIGHT_MODULE`, `CHROME`, and `FFMPEG` can override tool locations. GIF
export can also use the `imageio-ffmpeg` package. Rendering uses one browser
and four video-encoder threads. It makes no network requests and starts no
mapping searches. Preview mode exports selected review frames without a video.
