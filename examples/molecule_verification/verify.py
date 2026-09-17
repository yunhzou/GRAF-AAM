"""Verify inferred connectivity from two XYZ files using an unswept GRAFT match.

The bundled candidate is a controlled demonstration, not generative-model output.
XYZ has no bond labels: this example infers a binary graph independently on each
side with RDKit's covalent-radius rule. Bond order, charge, stereochemistry and
energetic stability are outside this connectivity check.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from graft import AAMProblem, AAMSearchConfig, MolecularEndpoint, search_aam
from graft.artifacts import write_aam_checkpoint


def endpoint_from_xyz(path):
    molecule = Chem.MolFromXYZBlock(Path(path).read_text())
    if molecule is None:
        raise ValueError(f'Cannot read XYZ: {path}')
    rdDetermineBonds.DetermineConnectivity(molecule, useVdw=True, covFactor=1.25)
    return MolecularEndpoint(tuple(a.GetSymbol() for a in molecule.GetAtoms()),
        molecule.GetConformer().GetPositions(), Chem.GetAdjacencyMatrix(molecule).astype(float),
        label=Path(path).stem)


def connection_check(problem, mapping):
    """Check both edges and nonedges under a complete element-preserving bijection."""
    n = problem.source_atom_count
    complete = (set(mapping) == set(range(n)) and
                sorted(mapping.values()) == list(range(problem.target_atom_count)) and
                all(problem.reactant.elements[i] == problem.product.elements[j] for i,j in mapping.items()))
    if not complete:
        return dict(complete=False, missing=None, extra=None, same_connectivity=False)
    ids = [mapping[i] for i in range(n)]
    r = problem.reactant.wbo >= .5
    p = problem.product.wbo[np.ix_(ids, ids)] >= .5
    missing = int(np.triu(r & ~p, 1).sum())
    extra = int(np.triu(p & ~r, 1).sum())
    return dict(complete=True, missing=missing, extra=extra, same_connectivity=not (missing or extra))


def verify(target, candidate, output, capture=False):
    from graft.growth.native import available
    target_endpoint, candidate_endpoint = endpoint_from_xyz(target), endpoint_from_xyz(candidate)
    problem = AAMProblem(target_endpoint, candidate_endpoint,
                         f'{target_endpoint.atom_count}-atom connectivity verification')
    config = AAMSearchConfig(sweep_cuts=False, seed_count=1, branch_limit=100,
                             iso_tolerance=.1, graph_floor=.5, cut_floor=.5)
    start = time.perf_counter()
    result = search_aam(problem, config, workers=1,
                        execution='reused_native' if available() else 'reference')
    elapsed = time.perf_counter()-start
    selected = None
    # A complete fragment alone is insufficient: reject extra candidate edges too.
    for terminal in result.graph.terminals:
        path = next(result.graph.paths(terminal))
        placements = [result.graph.transitions[e] for e in path.transitions
                      if result.graph.transitions[e].match is not None]
        if len(placements) != 1:
            continue
        mapping = dict(result.graph.states[terminal].mapping)
        check = connection_check(problem, mapping)
        if check['same_connectivity']:
            selected = (terminal, mapping, check)
            break
    output = Path(output);output.mkdir(parents=True, exist_ok=True)
    write_aam_checkpoint(result, output/'aam.checkpoint')
    report = dict(status='verified_connectivity' if selected else 'not_verified',
        atoms=problem.source_atom_count, target_edges=int(np.triu(problem.reactant.wbo,1).sum()),
        config=vars(config), search_seconds=elapsed, workers=1,
        input_sha256={Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (target,candidate)},
        connectivity_inference=dict(method='RDKit DetermineConnectivity', useVdw=True, covFactor=1.25),
        scope='Connectivity only; does not establish bond orders, stereochemistry, charge or stability.')
    if selected:
        terminal,mapping,check=selected
        report.update(terminal=terminal,fragments=1,mapping=mapping,**check)
        if capture:
            from graft.search_trajectory import build_trajectory
            trace=build_trajectory([dict(aam=result, context=0, terminals=[terminal])], title='Molecule verification: recorded atom-by-atom growth')
            (output/'trajectory.json.gz').write_bytes(gzip.compress(json.dumps(trace,separators=(',',':')).encode(),mtime=0))
    (output/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('mapping','input_sha256','config')},indent=2))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target',type=Path,default=Path(__file__).with_name('target.xyz'))
    parser.add_argument('--candidate',type=Path,default=Path(__file__).with_name('candidate.xyz'))
    parser.add_argument('--output',type=Path,default=Path('verification-output'))
    parser.add_argument('--capture',action='store_true')
    args=parser.parse_args()
    verify(args.target,args.candidate,args.output,args.capture)
