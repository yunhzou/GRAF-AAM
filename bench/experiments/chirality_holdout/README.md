# Chirality interpolation audit

See [the complete 140-case report](../../../reports/chirality_holdout_20260917/README.md)
for scope, reproduction commands, timings, and the offline viewer. The runner
loads trusted saved checkpoints; it does not rerun AAM or decoding. Use
`--minimum-only` to reproduce the original 166-candidate holdout viewer. Omit it
to inspect the full saved candidate list, with the same per-case watchdog.

For browser QA, set `PLAYWRIGHT_MODULE` to an installed Playwright module and,
if needed, `CHROME_PATH` to a Chrome executable, then run:

```sh
node bench/experiments/chirality_holdout/check_viewer.cjs /absolute/viewer.html /absolute/qa-output
```

Use a fresh output directory after changing the solver or source archives; the
runner's resume check distinguishes source directories and candidate scope, not
arbitrary edits within an existing source directory.
