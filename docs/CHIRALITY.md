# Chirality after signed-event decoding

This experimental API lives on `dev/chirality-postprocessing`; publication
`main` is unchanged. It reconnects coordinate/index-orientation selection to
the current saved AAM families. Search and bond-event decoding remain separate.

```python
from graft.postprocessing import decode_events
from graft.chirality import ChiralityConfig, select_chiral_witness

decoded = decode_events(aam)
candidate = decoded.candidates[0]
selected = select_chiral_witness(decoded, candidate)
if selected.status == "allowed":
    mapping = selected.mapping  # complete R -> P witness
    print(selected.diagnostics["ordinary_frames"])
else:
    print(selected.status, selected.diagnostics["reason"])
```

A self-contained executable example is in
[`examples/chirality/example.py`](../examples/chirality/example.py).

## Meaning of the constraint

This is **index orientation**, including permutations of chemically identical
ligands. It is not CIP assignment. The old orientation determinants and
permutation-invariant degeneracy rule are reused.

For each mapped center, find the reactant neighbors that remain bonded to its
product image at `graph_floor`. Three persistent ligands define a
center-relative orientation; four define an affine tetrahedron of the ligands.
Changes from four to five total ligands do not erase the persistent four-ligand
constraint. A planar or geometrically undefined frame carries no sign.

The default `mode="mutable"` preserves an ordinary frame when the saved AAM
union permits a nontrivial setwise ligand shuffle with its center fixed. This
mutability is proved using the actual correlated actions, not inferred from
independent atom orbits. An unavoidable fixed orientation change is reported as
`fixed_orientation_change`; it can represent a genuine reaction inversion.
For structure verification, use `mode="all"` to reject every defined ordinary
orientation reversal, including one at a center without assignment freedom.
Neither mode establishes E/Z stereochemistry or isotopic/CIP identity.

For higher coordination, ordinary hard constraints are installed first.
Defined source triples and four-ligand simplices are considered in decreasing
endpoint geometric robustness. Retain a frame only if it can coexist with all
previously retained frames and the ordinary constraints. Excluded frames are
reported in `reconfigured_high_coordinate_frames`. This is a **maximal feasible
basis**, not a maximum-cardinality optimum or evidence for a physical pathway.
`high_coordinate="strict"` instead requires all these frames. Lost connections
and undefined endpoint orientations make a frame inactive.

The old solver grouped dependent frames into automorphism orbits of one
analytical coset. The new saved relation can be an ordered product or a union
that is not a group. Consequently the new basis is defined over explicit local
source frames in that relation; it does not promise the old orbit grouping,
selected mapping. The ordinary constraint and its coordinate measurement are
retained. Both APIs select a feasible oriented witness without global geometry
ranking.

## What is certified

An allowed result is either the original saved witness, checked geometrically,
or a realization of an existing saved family's exact action program. It
preserves the candidate's **concrete signed broken/formed edges**, including
explicit hydrogens and the decoder's metal threshold. This keeps the selected
reaction core fixed; it is stricter than equality modulo the decoder's event
symmetry. It can intentionally exclude another member of the same event class
whose changed bonds involve different source atoms.

All saved families are eligible, including ones absent from the candidate's
support list because decoding was interrupted. No new symmetry is reconstructed
and no orbit is treated as an independent swap permission. Anchors encoded in
the saved families remain binding. The input AAM and decoded candidates are
not modified.

`forbidden` means no mapping in those saved complete families meets this policy,
not that the chemical reaction is impossible. `unknown` records an interrupted
proof and is never converted into a conflict. An original unchanged witness
has `family_id=None` and no new action certificate; a repaired witness carries
the realizing family ID and action sequence.

## Why this stays compact

The existing symbolic family compiler keeps ordered and coupled choices intact.
The selector checks a proposed witness and, when needed, adds a local orientation
constraint. One constraint covers all six orders of a three-ligand frame or all
24 orders of a four-ligand frame: target geometry is computed for the unordered
ligand set, and the parity of the ordered indices supplies its sign. The rule
is guarded by its actual center, ligand set and persistent connections.

Thus a failed frame excludes its entire wrong-parity local assignment class,
not just one full mapping. There is no list of all group elements or atom
bijections. Geometry, compiled families and mutability answers are cached for
the call. A witness that already satisfies all constraints can return without
compiling a solver. Solver checks and local refinements are reported explicitly.
Worst-case symbolic search is still combinatorial; these checks are not a
polynomial-time guarantee.

## Controls

| Setting | Default | Meaning |
|---|---|---|
| `graph_floor` | `0.2` | Persistent connectivity |
| `orientation_tolerance` | `0.1` | Ordinary/affine normalized-volume degeneracy |
| `group_orientation_tolerance` | `0.0` | High-coordinate triple degeneracy, with numerical error protection |
| `mode` | `"mutable"` | Historical shuffle-sensitive index orientation; `"all"` is stricter |
| `high_coordinate` | `"maximal"` | Report a maximal feasible basis; `"strict"` requires all frames |
| `seconds` | `None` | Optional soft time budget; no mapping-count or branch cap is added |

Use an external process watchdog for a hard limit. A timeout returns `unknown`
when observed inside the selector. Search limits remain the caller's existing
AAM configuration and are not changed here.

The selector returns a feasible mapping and orientation diagnostics, with no
RMSD calculation or ranking. The analytical `select_rp_mappings` API likewise
selects the first feasible saved branch/event coset, retaining its source
witness whenever feasible. Its fixed-correspondence rigid-fit diagnostic is
computed only after selection and cannot change the mapping.
