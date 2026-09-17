# Gold rearrangement — Grow. Branch. Decode.

A **22-second 3D film**, using the same presentation as the original Golden growth/branch/decode animation: matched molecules on the left, a growing tree on the right, one-second branch highlights, and two final candidate cards feeding the right-hand decoder.

![GRAFT: 3D fragment growth and two oxygen assignments in a gold-catalyzed rearrangement](gold-oxygen-preview.gif)

[Full-resolution MP4](gold-oxygen-3d.mp4) · [Offline interactive viewer](index.html) · [Inputs, witnesses and reproduction](../../../examples/gold_rearrangement/README.md)

Download and open `index.html` in a browser. It is self-contained and includes play/pause, scrubbing and “See candidates.”

## Follow one oxygen

The ring and **O*** label identify the original epoxide oxygen. Fragment colors show the atom correspondence throughout growth. The final cards enlarge the organic core and track that same O* into the product:

- **Candidate 1:** O* becomes the ester-link oxygen, consistent with literature route a.
- **Candidate 2:** O* becomes the ketone oxygen, consistent with literature routes b and c.

The paper favors route a energetically. Routes b and c converge at intermediate 14 and share this endpoint oxygen pattern. Recovering these endpoint correspondences supports different pathway hypotheses; GRAFT did not compute the intermediate mechanisms or their kinetics.

Source: González Pérez et al., “Mechanism of the Gold-Catalyzed Rearrangement of (3-Acyloxyprop-1-ynyl)oxiranes: A Dual Role of the Catalyst,” *J. Org. Chem.* **2009**, 74, 2982–2991. [DOI: 10.1021/jo802516k](https://pubs.acs.org/doi/10.1021/jo802516k).

## Timeline and evidence

**0–7 s:** the full 65-atom AuPPh₃ reactant/product complexes and two saved fragment-growth histories. **7–9 s:** blue and purple pulses highlight the two assignments of O*, with the camera paused. **9–14 s:** complete families feed into the decoder. **14–22 s:** two distinct candidate cards remain next to the decoder.

The tree shows two selected histories from different cut-sweep contexts, not the entire search. Every displayed node is a saved fragment commit; node numbers count assigned atoms. The final witnesses pass complete family-membership queries and have distinct canonical bond-event signatures. This is selected-witness decoding, not exhaustive decoding of all returned families. No bijection enumeration or fragment competition is used.

The final cards show the organic heavy-atom core. Red crosses mark breaking/weakening involving O*; green inward arrows mark forming/strengthening involving O*. Existing bonds remain black. The badge counts **all** signed threshold events in the full witness, including bond-order and metal changes; only the O*-centered changes are annotated visually. These counts are not elementary reaction steps.

Coordinates are supplied optimized endpoints, rigidly aligned for presentation. Camera rotations and core close-ups do not change molecular geometry. Dashed gold contacts display Au-pair WBOs between 0.15 and 0.45; this display rule is separate from matching and event thresholds. No endpoint morph or simulated dynamics is shown. The search settings, independently inferred oxygen provenance and GFN2-xTB input preparation are documented in the linked example.

## Rebuild

```bash
python manuscript/scripts/gold_film/build.py
node manuscript/scripts/gold_film/render.cjs manuscript/animations/gold_rearrangement 3d --video
python manuscript/scripts/gold_film/export_gif.py manuscript/animations/gold_rearrangement
```

The build rechecks family membership, event signatures and rigid display-coordinate transformations. Rendering requires Playwright, Chrome and ffmpeg; `PLAYWRIGHT_MODULE`, `CHROME`, and `FFMPEG` can override their locations. These commands reuse the frozen evidence, without rerunning AAM.

The MP4 is 1440 × 900 at 24 fps. The README GIF is 960 × 600 at 12 fps. `science-validation.json`, `rebuild-validation.json`, `browser-validation-3d.json` and `media-validation.json` record the checks. The earlier split 2D/3D presentation has been replaced; its history remains in Git.
