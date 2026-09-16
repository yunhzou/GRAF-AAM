# Final GRAFT end-to-end timing

Fresh attempts: 140/140; completed within the five-minute watchdog: 140.

One seed, single-edge cut sweep, cap 2,000, R-to-P, competition off. Event windows match the published 140-case comparison. No candidate-count cap. Four processes, each using one numerical thread.

## All attempted work, including timeouts

Observed process CPU spent: approximately 14.47 CPU minutes, or 6.20 CPU seconds per attempted reaction. This sampled process scope includes startup and interrupted calls. It is not the cost of completing all reactions.

**The following distributions cover only completed reactions. They are not full-dataset completed-runtime estimates when any reaction is interrupted.**

| Stage | Mean CPU s | Median CPU s | 95th percentile CPU s | Median wall s |
|---|---:|---:|---:|---:|
| Search | 1.257 | 0.300 | 3.468 | 0.330 |
| Final-family catalogue | 0.156 | 0.019 | 0.595 | 0.019 |
| Event decoding | 4.647 | 0.288 | 24.361 | 0.412 |
| Output serialization | 0.005 | 0.001 | 0.015 | 0.002 |
| End to end | 6.073 | 0.812 | 25.470 | 0.960 |

End to end includes input loading, search/checkpoint I/O, catalogue construction, decoding, candidate serialization, and light progress journaling. Initial imports and reference-set validation are excluded. The decoder includes lazy dependency initialization. Stage medians and percentiles are not additive.

Incomplete cases: []. See summary.json for resource-stop status, incurred process time, partial progress, and completed family certificates. An interrupted attempt is not a completed 300-second result.

All completed class sets and family counts match the archived baseline: True. No saved search or decoder outputs were reused in the timed pipeline.

Reproduction (from the repository root):

```sh
python native/build_engine.py
python bench/final_end_to_end.py --inputs /path/to/full140_inputs --output /path/to/new_run --cases all --workers 4
python bench/summarize_final_timing.py /path/to/new_run
```

The run manifest records source, native-binary and input hashes. Raw per-case search checkpoints, candidates and family journals remain in the run directory. This report preserves per-case timing and verification summaries.

The exact timed Python changes are preserved in `timed-source.patch`, against base commit `3749ded`. After timing, the empty-endpoint Boolean-mask dtype was made explicit and covered by a regression test; nonempty benchmark inputs already produce the identical Boolean arrays. The manifest retains the actual timed hashes.
