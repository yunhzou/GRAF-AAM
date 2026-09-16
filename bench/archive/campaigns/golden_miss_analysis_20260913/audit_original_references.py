"""Reparse original RDF labels independently of saved canonical atom indices."""
import json,hashlib,sys
from pathlib import Path
ROOT=Path(sys.argv[1]);ENGINE=Path(sys.argv[2]);OUT=Path(sys.argv[3]);sys.path[:0]=[str(ENGINE/'src'),str(ENGINE/'bench')]
import pynauty
from rdkit import Chem
from prepare_golden_benchmark import audit_block
from investigate_golden_mapping import original_reference_certificate
from golden_evaluation import colored_graph,project
rdf=ROOT/'golden-data/golden_dataset.rdf';blocks=rdf.read_text().split('$RFMT')[1:]
rows=[]
for case in [7,19,590,871,986,1228,1285,1358,1377,1475,1553]:
 original=audit_block(blocks[case],case);saved=json.loads((ROOT/f'golden-inputs/{case}/reference.json').read_text())
 verified=original_reference_certificate(original['mapped_reaction'])==pynauty.certificate(colored_graph(saved['features'],project(saved['mapping'],saved['features'])))
 assert verified
 sides=[]
 for smiles in original['mapped_reaction'].split('>>'):
  mol=Chem.MolFromSmiles(smiles)
  bonds={tuple(sorted((b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum()))):b.GetBondTypeAsDouble() for b in mol.GetBonds() if b.GetBeginAtom().GetAtomMapNum() and b.GetEndAtom().GetAtomMapNum()}
  sides.append(dict(bonds=bonds,labels={a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum()}))
 shared=sides[0]['labels']&sides[1]['labels']
 changes=[dict(labels=list(pair),weights=[s['bonds'].get(pair,0) for s in sides]) for pair in sorted(set(sides[0]['bonds'])|set(sides[1]['bonds'])) if set(pair)<=shared and sides[0]['bonds'].get(pair,0)!=sides[1]['bonds'].get(pair,0)]
 rows.append(dict(case=case,original_rdf_reference_verified=verified,mapped_reaction=original['mapped_reaction'],shared_labels=sorted(shared),changed_heavy_pairs=changes,reference_pair_count=len(shared)))
OUT.write_text(json.dumps(dict(rdf_sha256=hashlib.sha256(rdf.read_bytes()).hexdigest(),rows=rows),indent=2)+'\n')
print('Verified original RDF correspondence in',len(rows),'cases')
