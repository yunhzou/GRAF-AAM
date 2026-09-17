# Gold, oxygen and alternative atom histories

Two 30-second films of the **65-atom AuPPh₃-catalyzed rearrangement**, with independently recovered GRAFT witnesses. The 2D version emphasizes the organic core and oxygen identities. The 3D version shows the complete gold–triphenylphosphine complex using the supplied optimized coordinates.

![Two recovered oxygen-fate patterns in a gold-catalyzed rearrangement](gold-oxygen-preview.gif)

[2D film (MP4)](gold-oxygen-2d.mp4) · [3D film (MP4)](gold-oxygen-3d.mp4) · [Offline interactive viewer](index.html) · [Example, inputs and reproduction](../../../examples/gold_rearrangement/README.md)

Download `index.html` and open it in a browser. It is self-contained, works offline, and includes play/pause, scrubbing and a 2D/3D switch.

## Scientific reading

Oxygen A is the reactant epoxide oxygen, B the ester-link oxygen, and C the ester carbonyl oxygen. In recovered pattern A, oxygen A becomes the product ester-link oxygen and oxygen C the ketone oxygen. In pattern B, these two product assignments are reversed. Oxygen B becomes the ester carbonyl oxygen in both.

Pattern A is consistent with literature route a, starting with 1,3-ester migration. Pattern B is consistent with routes b and c, starting with 1,2-ester migration and epoxide activation. Those two routes converge at intermediate 14 and share the endpoint oxygen pattern. The study favors route a energetically. These are pathway hypotheses distinguished or constrained by endpoint correspondences; GRAFT did not calculate the intermediate mechanisms or their energetics.

Source: González Pérez et al., “Mechanism of the Gold-Catalyzed Rearrangement of (3-Acyloxyprop-1-ynyl)oxiranes: A Dual Role of the Catalyst,” *J. Org. Chem.* **2009**, 74, 2982–2991. [DOI: 10.1021/jo802516k](https://pubs.acs.org/doi/10.1021/jo802516k).

## What is animated

- **0–3 s:** gold-bound reactant 2 and product 9.
- **3–11 s:** six saved fragment commits for each of two selected histories from different sweep contexts. Node numbers count assigned atoms, including the catalyst and hydrogens omitted from the 2D drawing.
- **11–13 s:** the two selected mapping families blink in distinct colors.
- **13–20 s:** source oxygen identities are carried into the product; the 2D view adds red × for breaking/weakening and green inward arrows for forming/strengthening.
- **20–30 s:** two distinct oxygen-fate patterns side by side, linked to the literature hypotheses.

The tree is an excerpt, not the full 69-context search. The final two cards are certified, canonical-distinct witnesses selected from the deduplicated family catalogue; this film does not claim exhaustive event decoding. No atom bijection enumeration was used. Fragment competition is off.

2D bond orders are chemically assigned from the supplied organic geometries; the depiction suppresses hydrogens and stereochemical wedges for readability. Black lines retain existing bonds. 3D coordinates are rigidly aligned supplied endpoints, with a shared camera for final comparison. Dashed gold contacts display Au pair WBOs between 0.15 and 0.45; this display rule is separate from the matching and event thresholds. There is no morph between endpoints, simulated dynamics, or fabricated transition-state trajectory.

## Build and validate

```bash
python manuscript/scripts/gold_film/build.py
node manuscript/scripts/gold_film/render.cjs manuscript/animations/gold_rearrangement 2d --video
node manuscript/scripts/gold_film/render.cjs manuscript/animations/gold_rearrangement 3d --video
python manuscript/scripts/gold_film/export_gif.py manuscript/animations/gold_rearrangement
```

The build validates family membership and distinct event signatures. Rendering requires Playwright, Chrome and ffmpeg, with optional `PLAYWRIGHT_MODULE`, `CHROME`, and `FFMPEG` environment overrides. Only rendering is repeated by these commands; the frozen evidence is in `film-data.json`. Rerun the search separately using the example linked above.

`science-validation.json`, `rebuild-validation.json`, the two `browser-validation-*.json` files and `media-validation.json` document the checks. MP4s are 1600 × 900, 24 fps; the README GIF is 960 × 540, 12 fps. No source PDF or Gaussian output logs are bundled.
