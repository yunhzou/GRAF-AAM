# Canonical event decoding for the research film

The original film incorrectly treated three atom-indexed event lists as three unique event classes. The library's current `SignedEventIndex` assigns the same canonical ID to terminals 6085 and 6086: swapping source oxygen atoms O1 and O2 carries one signed pattern to the other. Terminal 6087 has a different canonical ID. The final film now groups equivalent event signatures and displays two decoded class representatives.

The fresh verification loads the saved one-seed, cap-2,000 sweep and competition archives for case 101, constructs `FinalBranchCatalogue`, rebuilds its shared symmetries, and symbolically decodes all 874 families from 411 unordered branches through five events. Both classes agree with the earlier complete saved result and with the fresh controlled cap comparison. No search or bijection enumeration is rerun.

Run `python replay_decode.py` from this directory in the repository Python environment to reproduce both class IDs from `final-catalogue.json.gz` and `input.json`. `final-decoding.json` records complete certificates for every family, representative mappings, fragment groups and source archive hashes. `witness-audit.json` shows the three original witnesses and the canonical ID collision. These are event classes on unverified mappings, not validated mechanisms.
