# GRAFT: Grow. Branch. Decode.

A 48-second offline 3D research preview. Open `index.html`, play `graft-grow-branch-decode.mp4`, or view the looping GIF. Black bonds show endpoint connectivity, pastel colors show matched fragments, red crosses mark bond loss/weakening, and green inward arrows mark formation/strengthening.

The recorded search tree has three terminal witnesses. A and B exchange equivalent oxygen atoms O1/O2 and have the same canonical event-pattern ID. The film groups them into one event class. C belongs to a second class. Different tree leaves are not automatically different final event candidates.

The final two cards use actual decoded witnesses from the deduplicated cap-2,000, one-seed sweep plus competition catalogue for coordinate case 101 (PR7): 411 fragment branches and 874 mapping families. A fresh symbolic decode is complete through five events and returns two canonical classes, each with four events. Source thresholds are 0.5 for ordinary pairs and 0.3 for metal-involving pairs. There is no explicit mapping-bijection or automorphism-group enumeration.

## Storyboard

- 0–7 s: grow a fragment with compatible placements.
- 7–19 s: retain two placements and grow under each condition.
- 19–28.5 s: show a second conditional fork and three complete tree leaves.
- 28.5–33 s: explain why A and B share an event class through O1/O2 exchange.
- 33–48 s: show two distinct event-class representatives from the complete final-catalogue decode.

## Evidence and scope

The growth illustration remains the actual archived context in `reports/pr7_search_trajectory_20260911/trace.json`, seed order 15 of a 30-order run, cap 1,000, isomorphism tolerance 1.0, and the C18–O21 sweep constraint. The finale uses the separate one-seed cap-2,000 catalogue for exactly the same endpoint WBOs, with its own provenance and complete event-window decode. It is not a claim that these three historical paths constitute that entire catalogue. Neither drawing nor mapping alternatives establish a physical reaction pathway. Coordinates undergo only rigid display transforms.

`reports/film_decoding_20260915/` provides the final catalogue, decoder results, every family certificate, the original endpoint input, and a replay script. `film-data.json` retains all three historical witnesses and their class membership, plus the two decoded final representatives. `science-validation.json` checks the canonical IDs and distinct final cards; browser validation checks timeline injection and final-card uniqueness.

## Rebuild

```sh
python manuscript/scripts/research_film/build_film.py --repo . --output manuscript/animations/graft_research_preview
node manuscript/scripts/research_film/render_film.cjs manuscript/animations/graft_research_preview video
python manuscript/scripts/research_film/export_gif.py manuscript/animations/graft_research_preview
```

The build uses the repository Python package, NumPy, Playwright, Chrome/Chromium and FFmpeg. `PLAYWRIGHT_MODULE`, `CHROME`, and `FFMPEG` can override paths. Rendering is offline and uses four encoder threads. `python reports/film_decoding_20260915/replay_decode.py` repeats the symbolic event-window decode without any new mapping search.
