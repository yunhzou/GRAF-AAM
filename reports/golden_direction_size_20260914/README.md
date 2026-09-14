# Golden recovery by endpoint size

Derived from final one- and two-seed saved directional verdicts, using all 1,851 input files verified against the original campaign SHA-256 manifest. No AAM searches or reference checks were rerun. Endpoint size counts explicit atoms including H. Equal-size cases use stored input order for smaller-first and reverse for larger-first. All directional verdicts at these settings are resolved. Per-case counts and input hashes, cohort summaries, and derivation sources are included in direction_recovery.json. The union is checked against the published bidirectional score.

Regenerate with: `python derive.py --repo /path/to/repo --inputs /path/to/golden-inputs --output direction_recovery.json`.
