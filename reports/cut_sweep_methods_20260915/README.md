# Cut-sweep Methods clarification

Added Section2.2, “Cut sweep to vary greedy fragment boundaries.” Verified against `src/rxn_core/aam.py::_search_cut` and `alignment/sweep.py::cut_sweep_items`: each subrun removes one source edge, preventing its use by growth and omitting its preservation test; input matrices remain intact. The uncut and single-edge results are combined. One seed means one ordering per subrun, including the entire sweep.

A single search perturbation can redirect downstream growth and expose a mapping with multiple bond events. The paper explicitly separates search-cut count from event count and does not suggest that discovering multi-event alternatives requires enumerating event combinations. Masks neither require the cut bond to break nor force its endpoints into different final fragments. The overall greedy search remains nonexhaustive.

Corrected the appendix's stale coordinate127 figure description to match the already-selected Golden case9, including reference provenance and H-event bounds. No algorithm changes or experiment reruns; benchmark results unchanged.

All21 pages rendered and visually reviewed. Numerical, molecular, reference and overflow validation passed.

PDF SHA256: f63bd5f66fadffffb8e28dfb66c166d9e8fe84312915571784825aa386c3598d
