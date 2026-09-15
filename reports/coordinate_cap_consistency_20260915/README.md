# Paired coordinate branch-cap study

See [RESULTS.md](RESULTS.md) for the completed 140-reaction cap-100 versus cap-2,000 comparison. Cap 100 is the default. Both runs use one forward seed ordering, the cut sweep, fragment competition, and separate event decoding under identical thresholds and fixed windows.

Cap 100 returns complete mappings for 138 reactions; cap 2,000 returns them for 140. Cases 123 and 125 are the only losses. The other 138 reactions have identical event-pattern sets across the whole decoded windows. Every cap-2,000 window reproduces its saved publication counterpart.

## Audit the saved evidence

Decompress `comparison.json.gz` to `comparison.json` in this directory. Run `python audit.py`, then `python report.py` to check the paired records and regenerate the report. `analysis.json` stores per-case differences and same-cohort stage timing. `saved-windows.json` freezes the earlier complete event-pattern sets; `controls.json` retains the prior SLAP comparison and fixed event windows.

## Reproduce searches

`manifest.json` identifies the exact source commit, source/native-source hashes, input hashes, and settings. Copy that revision's source tree into `engine/src` and compile its native extensions in the chosen Python environment. Supply the 140 matching endpoint records as `inputs/<case>/input.json`; the hashes in the manifest verify their identities. The input records contain `reactant` and `product` endpoint fields accepted by `MolecularEndpoint`. These are the same input records used by the current-validation coordinate study.

Set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `NUMEXPR_NUM_THREADS` to 1; use `PYTHONHASHSEED=0`. Run `python run.py case <case>` for each case from 0 through 139, then `python collect.py`. The driver saves unique attempt directories and refuses to overwrite a selected result. Its five-minute stage watchdog and memory guard preserve failed/partial work. Decoder continuation resumes completed family certificates; search is never resumed into a fresh-search timing.

`collect.py` adds internal search/competition CPU records to the basic aggregation and distinguishes empty searches from unfinished decoding. `aggregate.py` is retained unchanged from the frozen run; its older `unresolved` field meant that no pair of defined minima was available. The published `comparison.json.gz` uses the corrected classifications, and `analysis.json` is the final coverage summary. Search and competition timings exclude initial imports; worker timings include them. All decoder passes are included. Timings from different CPUs are kept separate.

Native binary identity was checked against the frozen native sources and verified in the three-case pilot; `pilot-proof.json` records its outcome and binary hashes. The small local smoke test is excluded from the measurements. These are unverified coordinate mappings: compare event counts and event-pattern coverage, not atom-mapping accuracy or mechanistic correctness.
