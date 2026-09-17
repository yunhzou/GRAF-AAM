"""Disconnect one bridge in the demonstration XYZ, checking the inferred graph."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds


def prepare(candidate, provenance, output):
    candidate=Path(candidate);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    m=Chem.MolFromXYZBlock(candidate.read_text());elements=[a.GetSymbol() for a in m.GetAtoms()]
    coords=m.GetConformer().GetPositions();order=json.loads(Path(provenance).read_text())['permutation_new_to_original']
    # Original target bond C33-C37, with zero-based indices 32 and 36.
    a,b=order.index(32),order.index(36)
    def text(x):
        return f'{len(x)}\nControlled negative example: one disconnected C-C bridge; not a reaction trajectory\n'+''.join(f'{e} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n' for e,p in zip(elements,x))
    def adjacency(x):
        mol=Chem.MolFromXYZBlock(text(x));rdDetermineBonds.DetermineConnectivity(mol,useVdw=True,covFactor=1.25)
        return Chem.GetAdjacencyMatrix(mol)
    original=adjacency(coords);g=nx.from_numpy_array(original)
    assert nx.is_connected(g) and frozenset((a,b)) in {frozenset(e) for e in nx.bridges(g)}
    g.remove_edge(a,b);moving=sorted(nx.node_connected_component(g,b))
    assert a not in moving
    direction=coords[b]-coords[a];direction/=np.linalg.norm(direction)
    if len(moving)>len(coords)/2:
        moving=sorted(set(g)-set(moving));direction=-direction
    # Choose the smallest tested separation that removes only this bridge.
    for displacement in np.arange(3.,21.,.5):
        trial=coords.copy();trial[moving]+=displacement*direction;actual=adjacency(trial)
        expected=original.copy();expected[a,b]=expected[b,a]=0
        if np.array_equal(actual,expected):break
    else:raise RuntimeError('Could not isolate exactly one broken connection')
    (output/'candidate-broken.xyz').write_text(text(trial))
    record=dict(kind='controlled_disconnection',input_sha256=hashlib.sha256(candidate.read_bytes()).hexdigest(),
        source_target_edge=[32,36],candidate_edge=[a,b],moving_atoms=moving,
        translation=(displacement*direction).tolist(),translation_distance=float(displacement),
        original_edges=int(original.sum()/2),broken_edges=int(actual.sum()/2),
        original_components=1,broken_components=nx.number_connected_components(nx.from_numpy_array(actual)),
        removed_connections=[[a,b]],added_connections=[],
        scope='Coordinate perturbation for verification failure; not a simulated bond-breaking mechanism.')
    (output/'break-preparation.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,default=Path(__file__).with_name('candidate.xyz'))
    p.add_argument('--provenance',type=Path,default=Path(__file__).with_name('preparation.json'))
    p.add_argument('--output',type=Path,default=Path(__file__).parent)
    a=p.parse_args();prepare(a.candidate,a.provenance,a.output)
