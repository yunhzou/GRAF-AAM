"""Regenerate the controlled demo candidate from the bundled target XYZ."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from rdkit import Chem,rdBase
from rdkit.Chem import rdDetermineBonds
import networkx as nx


def prepare(target, output):
    target=Path(target);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    molecule=Chem.MolFromXYZBlock(target.read_text())
    elements=[a.GetSymbol() for a in molecule.GetAtoms()]
    xyz=molecule.GetConformer().GetPositions();n=len(elements)
    def xyz_text(c,e=elements,title=''):
        return str(n)+'\n'+title+'\n'+''.join(f'{s} {x:.8f} {y:.8f} {z:.8f}\n' for s,(x,y,z) in zip(e,c))
    def adjacency(c):
        m=Chem.MolFromXYZBlock(xyz_text(c))
        rdDetermineBonds.DetermineConnectivity(m,useVdw=True,covFactor=1.25)
        return Chem.GetAdjacencyMatrix(m)
    def rotate(c,origin,axis,angle):
        axis=axis/np.linalg.norm(axis);p=c-origin;t=np.deg2rad(angle)
        return origin+p*np.cos(t)+np.cross(axis,p)*np.sin(t)+np.outer(p@axis,axis)*(1-np.cos(t))
    original=adjacency(xyz);graph=nx.from_numpy_array(original)
    candidate=xyz.copy();changes=[]
    # These bridge indices belong to the bundled holdout-68 target, not a
    # general conformer generator. Changes are retained only if connectivity stays fixed.
    for u,v,angle in [(2,10,30),(32,36,-22),(33,37,26),(34,38,-24),(35,39,28)]:
        cut=graph.copy();cut.remove_edge(u,v);part=sorted(nx.node_connected_component(cut,v))
        if len(part)>n/2:part=sorted(set(graph)-set(part))
        trial=candidate.copy();trial[part]=rotate(candidate[part],candidate[u],candidate[v]-candidate[u],angle)
        if np.array_equal(adjacency(trial),original):
            candidate=trial;changes.append(dict(bond=[u,v],degrees=angle,atoms=part))
    rng=np.random.default_rng(20260917);noise=rng.normal(0,.012,candidate.shape)
    noise_applied=np.array_equal(adjacency(candidate+noise),original)
    if noise_applied:candidate+=noise
    order=rng.permutation(n)
    candidate=rotate(candidate,np.zeros(3),np.array([.3,.7,1.]),137)+[4,-3,2]
    assert np.array_equal(adjacency(candidate),original)
    (output/'candidate.xyz').write_text(xyz_text(candidate[order],[elements[i] for i in order],
        'Controlled demonstration candidate; conformational perturbations, rigid rotation and shuffled atom order; not model output'))
    (output/'preparation.json').write_text(json.dumps(dict(source_case=68,source_name='pr17.carbene.ins_ts6a',
        endpoint='reactant',target_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        atom_count=n,heavy_atoms=sum(e!='H' for e in elements),
        candidate='Controlled demonstration, not output from a generative model',random_seed=20260917,
        torsion_changes=changes,coordinate_noise_sigma=.012,noise_applied=bool(noise_applied),
        permutation_new_to_original=order.tolist(),rdkit_version=rdBase.rdkitVersion,
        inference=dict(method='RDKit DetermineConnectivity',useVdw=True,covFactor=1.25),
        target_edges=int(original.sum()/2)),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target',type=Path,default=Path(__file__).with_name('target.xyz'))
    parser.add_argument('--output',type=Path,default=Path(__file__).parent)
    args=parser.parse_args();prepare(args.target,args.output)
