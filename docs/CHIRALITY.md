# Chirality after signed-event decoding

Experimental API on `dev/chirality-postprocessing`. Search and event decoding remain separate. The default is a **strict coordinate-orientation filter over all saved AAM families**, followed by selection of one feasible witness. It does not optimize global RMSD, interpolate coordinates, or rank collision-free paths.

```python
from graft.chirality import select_chiral_witness, query_chiral_witness

candidate = decoded.minimum_candidates[0]
result = select_chiral_witness(decoded, candidate)
if result.status == "allowed":
    mapping = result.mapping
elif result.status == "forbidden":
    print("No orientation-preserving mapping exists in these saved AAM families")
else:  # unknown: watchdog/solver could not establish the answer
    print(result.diagnostics["reason"])

# Test an entire old mapping, or ask whether a few assignments can coexist.
check = query_chiral_witness(decoded, candidate, old_mapping)
conditional = query_chiral_witness(decoded, candidate, {source_atom: product_atom})
```

The query is an implicit representation of the filtered subset. No bijections or group closures are expanded. Selecting one witness does not discard other feasible families. A preferred witness can be retained with `select_chiral_witness(decoded, candidate, preferred_mapping=old_mapping)`; it is used only after certification of membership, exact events, and orientation. A preference is not a hard constraint; use the query to require particular assignments.

## What the filter means

At bond floor 0.2, find the reactant neighbors that remain bonded to each mapped product center. Three persistent neighbors define a center-relative signed volume. Four define the signed volume of their affine tetrahedron. Require sign preservation when both endpoints are geometrically defined.

At centers with more than four neighbors **at either actual mapped endpoint**, also require all defined persistent three- and four-neighbor frame orientations to agree. An unrelated high-coordinate atom of the same element does not activate extra constraints on a four-coordinate center. Four-neighbor frames use a permutation-invariant volume normalization; ordinary degeneracy tolerance is 0.1, and the high-coordinate triple tolerance is 0.0 with a numerical roundoff guard. These tolerances are configurable.

A frame with a lost connection or undefined endpoint orientation is inactive. The filter does not assign CIP or E/Z labels, distinguish isotopes absent from the input, or establish that a proposed reaction must preserve stereochemistry. Real reactions can invert centers or rearrange coordination. Consequently, a strict conflict is a statement about these orientation requirements and the **saved** AAM search results, not proof that the reaction or every possible atom mapping is impossible.

All decoded candidates keep their concrete broken and formed edges, including the selected metal threshold. This is stronger than preserving an event class up to endpoint symmetry. Saved anchors and correlated, ordered symmetry actions remain authoritative.

## Results and infeasibility

- `allowed`: a complete saved-family witness satisfies every applicable requested orientation constraint. Membership is certified even if the preferred witness already looks geometrically valid.
- `forbidden`: every saved family was ruled out for the requested assignments, events, and orientation constraints. `search_exhaustive=True` records this scope. There is no replacement mapping.
- `unknown`: a timeout or inconclusive solver result; never interpreted as an empty subset.
- `relaxed`: an explicitly requested historical relaxation returned a mapping with full-orientation violations. Its violations are listed; it is **not** a strict successful correction.

An `allowed` result is not a clash-free interpolation guarantee. Even several mappings that preserve all local signs can differ in torsion, fragment placement, and interpolation quality. Those are separate questions. The selected witness is first feasible, with an optional user preference; it is not claimed to be the best geometric representative.

The default selector and subset query use the same stable predicate, independent of whichever witness happens to be found first. In particular, they do not greedily choose one high-coordinate frame basis, discard competing bases, and then present that choice as the entire allowed subset.

## Explicit historical relaxation

```python
from graft.chirality import ChiralityConfig

historical = ChiralityConfig(
    mode="mutable",
    high_coordinate="maximal",
    high_coordinate_scope="selected_family",
)
relaxed = select_chiral_witness(decoded, candidate, historical)
print(relaxed.status, relaxed.diagnostics.get("orientation_violations", []))
```

`mutable` permits fixed ordinary orientation reversals when the saved relation has no setwise ligand shuffle. `maximal` greedily retains a feasible high-coordinate frame basis and reports omitted frames. It is a heuristic relaxation, not a complete representation of all admissible reconfiguration choices. It remains opt-in for comparison with older work; subset queries reject these reference-dependent policies. The old viewer's “zero selected violations” excluded deliberately relaxed frames, so it was not evidence of full preservation.

## Configuration

| Option | Default | Meaning |
|---|---|---|
| `mode` | `"all"` | All defined ordinary orientations; `"mutable"` explicitly relaxes fixed changes |
| `high_coordinate` | `"strict"` | All applicable higher-coordinate frames; `"maximal"` is historical relaxation |
| `high_coordinate_scope` | `"union"` | Only controls the explicit maximal relaxation; strict always searches the full union |
| `graph_floor` | `0.2` | Endpoint connectivity threshold |
| `orientation_tolerance` | `0.1` | Ordinary and four-neighbor frame degeneracy threshold |
| `group_orientation_tolerance` | `0.0` | Higher-coordinate three-neighbor frame threshold |
| `seconds` | `None` | Optional soft watchdog; use isolated processes for a hard limit |

Rejection bounds use singleton images of the ordered family action program. They can prove an unavoidable orientation conflict but never certify a mapping. All positive answers use the full correlated solver. Exhausted family models are released to keep memory bounded by active work; this does not cap the number of families searched.

## Validation

Tests cover exact and partial subset membership, cross-family alternatives, unchanged inputs, fixed inversions, conflicting correlated permutations, exact event preservation, anchors, proper rotations and relabeling, planar frames, solver timeouts, explicit relaxed statuses, and finite exhaustive oracles for small symmetry programs. The example in `examples/chirality/example.py` is self-contained.

The old/new holdout audit is recorded in `reports/chirality_strict_audit_20260917/`. Earlier `reports/chirality_holdout_20260917/` results used the historical mutable/maximal policy and are superseded for strict chirality claims.
