# GRAFT: Grow. Compete. Decode.

A 28-second offline 3D preview using **Golden case 15**. It shows fragment growth, a conditional fork, an actual fragment-competition completion, and animated merging into canonical event classes on the right.

The original branches stay visible while the O18 fragment takes priority across C16. The blue return arrow opens a new completion with 28 retained assignments; O17 and H19 are refilled. At decoding, B and C flow into the same event class. Final cards show two distinct decoded representatives, including a Golden reference-equivalent class. The saved catalogue contains nine classes through six events; two are illustrated in detail.

Black bonds preserve endpoint connectivity; pastel colors identify matched fragments. Red crosses mark breaking/weakening and green inward arrows mark formation/strengthening. Coordinates are illustrative conformers, not reaction dynamics.

## Timeline

- 0–7 s: grow two fragments under successive conditions.
- 7–9 s: pulse branch 1 in blue for one second, then branch 2 in purple for one second. Changed atom correspondences, tree edges and nodes highlight together while the camera stays fixed.
- 9–11 s: complete the retained baseline branches A and B.
- 11–16.8 s: animate fragment competition, release conflicting assignments and complete C.
- 16.8–19.2 s: transition the right-hand tree into the decoder; B and C visibly merge.
- 19.2–28 s: show decoded molecule pairs and event lists, while retaining the merge diagram.

The [Golden evidence](../../../reports/golden_film_20260915/README.md) contains the portable archives, exact competition proposal, complete catalogue decode and a verification script. The growth tree is a selected excerpt, not the entire sweep. Competition is a new completion under revised priorities; the original alternatives are retained.

## Rebuild

From the repository root:

```sh
python manuscript/scripts/research_film/build_film.py --repo . --output manuscript/animations/graft_research_preview
node manuscript/scripts/research_film/render_film.cjs manuscript/animations/graft_research_preview video
python manuscript/scripts/research_film/export_gif.py manuscript/animations/graft_research_preview
```

The build uses NumPy, the repository package, Playwright, Chrome/Chromium and FFmpeg. `PLAYWRIGHT_MODULE`, `CHROME`, and `FFMPEG` can override executable/module paths. Rendering is offline and uses four encoder threads. Source mappings and provenance are in `film-data.json`; numerical, browser and media checks accompany the files.
