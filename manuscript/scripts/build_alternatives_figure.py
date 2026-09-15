"""Vector figure of current GRAFT candidates, checked from a six-atom API run."""
import json,os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).resolve().parents[1]/"build/mpl-cache"))
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from build_molecule_figure import molecule,INK,MUTED,TEAL,ORANGE,PURPLE,RED,BLUE,LINE

R='[C:1]([O:2][H:6])([H:3])([H:4])[H:5]'
A='[C:1](=[O:2])([H:3])[H:4]'
B='[C:1](=[O:2])([H:3])[H:6]'

def check_example(man):
 d=json.loads((man/'evidence/multicandidate_example.json').read_text())
 assert d['complete_saved_families'] and [c['total'] for c in d['candidates']]==[4,6]
 assert d['queries']=={'same_class_H3_H4_swap':'allowed','other_class_preserve_events':'forbidden','other_class_any_saved_events':'allowed'}
 wr=np.array(d['wbo_R']);wp=np.array(d['wbo_P'])
 expected=[[(0,1,1),(0,4,-1),(1,5,-1),(4,5,1)],[(0,1,1),(0,3,-1),(0,4,-1),(0,5,1),(1,5,-1),(3,4,1)]]
 from rdkit import Chem
 def bonds(s):
  p=Chem.SmilesParserParams();p.removeHs=False;m=Chem.MolFromSmiles(s,p);assert m
  return {tuple(sorted((b.GetBeginAtom().GetAtomMapNum()-1,b.GetEndAtom().GetAtomMapNum()-1))):b.GetBondTypeAsDouble() for b in m.GetBonds()}
 r=bonds(R)
 assert len(r)==5 and all(wr[i,j]==v for (i,j),v in r.items())
 for m,expected_events,smiles in zip(d['figure_witnesses_zero_based'],expected,[A+'.[H:5][H:6]',B+'.[H:4][H:5]']):
  m={int(k):v for k,v in m.items()}
  delta=wp[np.ix_([m[i] for i in range(6)],[m[i] for i in range(6)])]-wr
  events=[(i,j,int(np.sign(delta[i,j]))) for i in range(6) for j in range(i+1,6) if abs(delta[i,j])>=.5]
  assert events==expected_events
  b=bonds(smiles)
  assert len(b)==4
  assert all(b.get((i,j),0)==wp[m[i],m[j]] for i in range(6) for j in range(i+1,6))
 counts=json.loads((man/'evidence/output_multiplicity.json').read_text())
 assert all(m['cases']==1851 for m in counts['methods'])
 assert all(m['multiple_record_cases']==0 for m in counts['methods'] if not m['method'].startswith('slap'))
 assert all(m['multiple_record_cases']>0 for m in counts['methods'] if m['method'].startswith('slap'))
 return d


def check_case(man):
 from rdkit import Chem
 d=json.loads((man/'evidence/coordinate_case127.json').read_text())
 wr=np.array(d['input']['reactant']['wbo']);wp=np.array(d['input']['product']['wbo'])
 assert len(d['comparator_minimum_ids']['native_slap'])==1
 assert d['case']==127 and d['minimum']==4 and d['saved_window_complete']
 assert len(d['patterns'])==3 and len({p['id'] for p in d['patterns']})==3
 assert {p['id'] for p in d['patterns']}==set(d['comparator_minimum_ids']['slap_sweep'])
 for p in d['patterns']:
  v=p['mapping'];delta=wp[np.ix_(v,v)]-wr
  events={kind:[[i,j] for i in range(23) for j in range(i+1,23) if sign*delta[i,j]>=.5] for kind,sign in [('broken',-1),('formed',1)]}
  assert events==p['events'] and sum(map(len,events.values()))==p['total']==4
  assert sorted(a for f in p['source_fragments'] for a in f)==list(range(23))
  assert p['query_status']=='recovered'
  floor,iso=p['supporting_family']['policy']
  assert all(wp[v[i],v[j]]>=floor and abs(wp[v[i],v[j]]-wr[i,j])<=iso+1e-9 for i,j in p['supporting_family']['required_edges'])
 for side,w in [('reactant',wr),('product',wp)]:
  params=Chem.SmilesParserParams();params.removeHs=False;params.sanitize=False
  molecule=Chem.MolFromSmiles(d[side+'_smiles'],params);assert molecule.GetNumAtoms()==23
  edges={tuple(sorted((b.GetBeginAtom().GetAtomMapNum()-1,b.GetEndAtom().GetAtomMapNum()-1))):b.GetBondTypeAsDouble() for b in molecule.GetBonds()}
  assert edges=={(i,j):(2 if w[i,j]>=1.4 else 1) for i in range(23) for j in range(i+1,23) if w[i,j]>=.5}
 assert d['shuffle']['status']=='allowed' and d['shuffle']['same_event_id']==d['patterns'][0]['id']
 return d

# Same endpoint geometry for all rows; product identities change with the map.
COORDS={1:(0,1.5),6:(0,3),7:(-1.3,.75),4:(-1.3,-.75),3:(0,-1.5),2:(1.3,-.75),8:(1.3,.75),13:(-2.6,1.5),20:(-3.9,.75),9:(-2.6,-1.5),16:(2.6,-1.5)}
PALETTE=[(.70,.87,.81),(.98,.81,.49),(.81,.68,.92)]

def build(man):
 from rdkit import Chem
 from matplotlib.patches import Circle
 d=check_case(man);check_example(man)
 def smiles(s,mapping=None):
  params=Chem.SmilesParserParams();params.removeHs=False;params.sanitize=False
  mol=Chem.MolFromSmiles(s,params)
  Chem.SanitizeMol(mol,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
  if mapping:
   inv={v+1:k+1 for k,v in enumerate(mapping)}
   for atom in mol.GetAtoms():atom.SetAtomMapNum(inv[atom.GetAtomMapNum()])
  mol=Chem.RemoveHs(mol,sanitize=False)
  return Chem.MolToSmiles(mol,kekuleSmiles=True)
 fig=plt.figure(figsize=(10.8,9.23));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1296),ylim=(1108,0));ax.axis('off')
 def t(x,y,s,size=11,color=INK,weight='normal',ha='left'):
  ax.text(x,y,s,fontsize=size,color=color,weight=weight,ha=ha,va='center',linespacing=1.35,zorder=9)
 def box(x,y,w,h,fc='#F5F7FB',ec=None):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=12',facecolor=fc,edgecolor=ec or fc,lw=.8,zorder=0))
 def arrow(a,b,color=INK,rad=0,lw=1.5):
  ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,color=color,lw=lw,connectionstyle=f'arc3,rad={rad}',zorder=6))
 def mol(s,x,y,w,h,**kw):
  return molecule(ax,s,x,y,w,h,font=29,explicit_hydrogens=True,align_acetyl=False,preserve_bond_orders=True,**kw)
 def heading(x,y,letter,title):t(x,y,letter,15,weight='bold');t(x+31,y,title,13,weight='bold')
 box(12,12,1272,754)
 heading(32,40,'a','Case 127: three different mappings, the same minimum event count')
 t(32,70,'140-reaction coordinate set  ·  GRAFT: 3 classes  ·  SLAP sweep: 3  ·  Native SLAP: 1',10.5,MUTED)
 for x,color in zip([45,70,95],PALETTE):ax.add_patch(Circle((x,107),10,facecolor=color,edgecolor='none'))
 t(118,107,'Matched fragments keep their colors',10,MUTED)
 ax.plot([603,666],[107,107],color=RED,lw=2.5,linestyle=(0,(2.6,1.6)));t(680,107,'Break / weaken',10,RED)
 ax.plot([949,1012],[107,107],color=BLUE,lw=2.5);t(1026,107,'Form / strengthen',10,BLUE)
 for index,p in enumerate(d['patterns']):
  y=140+202*index;box(28,y,1238,194,'white','#DCE3E9')
  letter='ABC'[index];t(48,y+27,letter,15,TEAL,'bold');t(93,y+27,'4 events',11,TEAL,'bold')
  owners={a+1:PALETTE[k] for k,f in enumerate(p['source_fragments']) for a in f}
  changes=lambda key:{tuple(a+1 for a in pair):(-1 if key=='broken' else 1) for pair in p['events'][key]}
  # Fixed six-member perimeter reveals the extra reactant bridge clearly.
  mol(smiles(d['reactant_smiles']),183,y+9,425,179,owners=owners,bond_changes=changes('broken'),coordinates=COORDS)
  arrow((643,y+101),(718,y+101),INK,lw=1.8)
  inv={v+1:k+1 for k,v in enumerate(p['mapping'])}
  target_coords={inv[k]:xy for k,xy in COORDS.items()}
  mol(smiles(d['product_smiles'],p['mapping']),753,y+9,425,179,owners=owners,bond_changes=changes('formed'),coordinates=target_coords)
  t(96,y+91,'2 negative\n2 positive',9.1,MUTED,ha='center')
 t(647,751,'Atom numbers track reactant identities. Colors are assigned within each row; spectator H atoms are implicit.',9.8,MUTED,ha='center')
 box(12,783,624,311)
 heading(32,811,'b','Symmetry remains inside one candidate')
 t(32,843,'Candidate A: exchange two H atoms of methyl C₉.',10.3,MUTED)
 methyl='[*:24][C:9]([H:10])([H:11])[H:12]'
 owners={i:PALETTE[1] for i in [9,10,11,12]}
 left=mol(methyl,37,865,215,145,owners=owners,notes={i:i for i in [9,10,11,12]})
 right=mol(methyl,381,865,215,145,owners=owners,notes={9:9,10:11,11:10,12:12})
 # Trace individual H identities while preserving the fragment background color.
 for source,dest,color,rad in [(10,11,PURPLE,-.12),(11,10,BLUE,.12)]:
  a,b=left[source],right[dest]
  for center in [a,b]:ax.add_patch(Circle(center,17,facecolor='none',edgecolor=color,lw=1.5,zorder=8))
  arrow((a[0]+20,a[1]),(b[0]-20,b[1]),color,rad,lw=1.4)
 t(321,1041,'Certified allowed · same event class',11,TEAL,'bold','center')
 t(321,1070,'R = the rest of the molecule; all other assignments fixed.',9.2,MUTED,ha='center')
 box(650,783,634,311)
 heading(670,811,'c','Retain the fragment for geometric use')
 t(670,843,'Candidate B: the ethyl group stays together.',10.3,MUTED)
 ethyl='[*:24][CH2:13][CH3:20]';owners={13:PALETTE[2],20:PALETTE[2]}
 mol(ethyl,675,882,214,123,owners=owners,notes={13:13,20:20})
 mol(ethyl,1051,873,214,140,owners=owners,notes={13:13,20:20},rotation=-45)
 arrow((907,945),(1033,945),PURPLE,rad=-.35,lw=1.8)
 t(966,889,'rotate + translate',10,PURPLE,'bold','center')
 t(966,1041,'Preserve internal bonds; refine reacting contacts',10.5,TEAL,'bold','center')
 t(966,1070,'Possible TS preparation; no pathway validation is claimed.',9.2,MUTED,ha='center')
 for ext in ['pdf','svg','png']:fig.savefig(man/'figs'/f'fig4_alternatives.{ext}',facecolor='white',bbox_inches='tight',pad_inches=.025,dpi=240)
 plt.close(fig)

if __name__=='__main__':build(Path(__file__).resolve().parents[1])
