"""Self-contained coordinate/index-chirality selection, using only public APIs.

Run with: python examples/chirality/example.py
Requires the package's postprocessing extra. No external molecule or xTB needed.
"""
import numpy as np
from graft import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam
from graft.postprocessing import decode_events
from graft.chirality import ChiralityConfig, select_chiral_witness

# CH3F: three identical H ligands make index-orientation shuffles possible.
elements = ('C', 'F', 'H', 'H', 'H')
reactant_xyz = np.array([
    [0, 0, 0], [-.9, -.9, -.9], [.9, .9, -.9],
    [.9, -.9, .9], [-.9, .9, .9],
])
product_xyz = reactant_xyz.copy()
product_xyz[[3, 4]] = product_xyz[[4, 3]]
wbo = np.zeros((5, 5))
wbo[0, 1:] = wbo[1:, 0] = 1.0
problem = AAMProblem(
    MolecularEndpoint(elements, reactant_xyz, wbo),
    MolecularEndpoint(elements, product_xyz, wbo),
)
aam = search_aam(problem, AAMSearchConfig(sweep_cuts=False, branch_limit=200), workers=1)
decoded = decode_events(aam)
for candidate in decoded.minimum_candidates:
    selected = select_chiral_witness(decoded, candidate, ChiralityConfig(mode='all'))
    assert selected.status == 'allowed'
    assert decoded.query(candidate, selected.mapping).status == 'allowed'
    print('Bond events:', candidate.total)
    print('Selected mapping:', selected.mapping)
    print('Orientation frames:', selected.diagnostics['ordinary_frames'])
