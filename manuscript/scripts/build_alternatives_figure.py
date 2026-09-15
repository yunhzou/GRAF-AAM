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

def build(man):
 d=check_example(man)
 fig=plt.figure(figsize=(10.8,7.25));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1296),ylim=(870,0));ax.axis('off')
 def t(x,y,s,size=11,color=INK,weight='normal',ha='left'):
  ax.text(x,y,s,fontsize=size,color=color,weight=weight,ha=ha,va='center',linespacing=1.4,zorder=7)
 def box(x,y,w,h,fc='#F5F7FB',ec=None):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=13',facecolor=fc,edgecolor=ec or fc,lw=.9,zorder=0))
 def arrow(a,b,color=INK,rad=0):
  ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,color=color,lw=1.5,connectionstyle=f'arc3,rad={rad}',zorder=1))
 def mol(s,x,y,w,h,**kw):
  return molecule(ax,s,x,y,w,h,font=30,explicit_hydrogens=True,align_acetyl=False,**kw)
 def heading(x,y,letter,title):
  t(x,y,letter,15,weight='bold');t(x+32,y,title,13,weight='bold')
 box(12,12,1272,472)
 heading(32,43,'a','One endpoint pair, two GRAFT bond-event hypotheses')
 t(32,75,'Methanol → formaldehyde + H₂     |     All H atoms explicit; numbers track reactant identities',10.5,MUTED)
 mol(R,25,146,278,244)
 t(164,421,'Same starting atoms',11,weight='bold',ha='center')
 arrow((305,267),(350,267),TEAL)
 box(354,221,190,92,'#E3F0F0','#B1D5D4')
 t(449,245,'GRAFT',13,TEAL,'bold','center');t(449,274,'2 event classes',11,TEAL,'bold','center')
 t(449,299,'1 seed · sweep (default)',8.7,MUTED,ha='center')
 ax.plot([545,578,578],[267,267,180],color=TEAL,lw=1.5)
 ax.plot([578,578],[267,379],color=TEAL,lw=1.5)
 arrow((578,180),(620,180),TEAL);arrow((578,379),(620,379),ORANGE)
 for y,color,label,smiles,h2,lines in [
  (102,TEAL,'A  ·  4 changes',A,'[H:5][H:6]',('Break C₁–H₅, O₂–H₆','Form H₅–H₆; strengthen C₁–O₂')),
  (291,ORANGE,'B  ·  6 changes',B,'[H:4][H:5]',('Break C₁–H₄, C₁–H₅, O₂–H₆','Form C₁–H₆, H₄–H₅; strengthen C₁–O₂'))]:
  box(627,y,637,174,'white','#D9E1E8')
  t(646,y+24,label,11.5,color,'bold')
  mol(smiles,650,y+46,156,117)
  t(819,y+102,'+',17,MUTED,ha='center')
  mol(h2,833,y+66,96,62)
  t(947,y+68,lines[0],8.4,RED)
  # Break second long line to keep text inside the outcome cards.
  words=lines[1].replace('; strengthen',';\nstrengthen')
  t(947,y+111,words,8.4,BLUE)
 box(12,501,624,354)
 heading(32,533,'b','A representative still has allowed shuffles')
 t(32,565,'Class A: exchange the two retained methyl H atoms',10.3,MUTED)
 mol(A,42,591,182,149)
 arrow((241,660),(350,660),PURPLE)
 t(295,620,'H₃ ↔ H₄',12,PURPLE,'bold','center')
 mol(A,365,591,182,149,notes={1:1,2:2,3:4,4:3})
 t(295,742,'H₅–H₆ is unchanged',10.3,MUTED,ha='center')
 t(34,781,'Same-class query: allowed',11,TEAL,'bold')
 t(34,817,'A → B: forbidden in A; allowed across saved classes',9.8,MUTED)
 box(650,501,634,354)
 heading(670,533,'c','Fragments can guide geometric preparation')
 t(670,565,'Illustrative retained fragment: move C–C–O together',10.3,MUTED)
 # The same connected motif, depicted in two orientations. This panel is a
 # downstream concept, not a fragment extracted from the six-atom example.
 atoms={7:(.73,.88,.87),8:(.73,.88,.87),9:(.73,.88,.87)}
 mol('[CH3:7][CH2:8][OH:9]',669,611,206,112,owners=atoms,notes=False)
 mol('[CH3:7][CH2:8][OH:9]',1054,597,206,140,owners=atoms,notes=False,rotation=-45)
 arrow((886,668),(1040,668),TEAL,rad=-.3)
 t(963,611,'rotation + translation',9.5,TEAL,'bold','center')
 t(966,744,'Keep internal bonds; refine reactive contacts',10,INK,ha='center')
 t(670,781,'Then test a transition-state guess',11,TEAL,'bold')
 t(670,817,'Proposed use · no barrier or pathway validation here',9.8,MUTED)
 for ext in ['pdf','svg','png']:
  fig.savefig(man/'figs'/f'fig4_alternatives.{ext}',facecolor='white',bbox_inches='tight',pad_inches=.025,dpi=240)
 plt.close(fig)

if __name__=='__main__':build(Path(__file__).resolve().parents[1])
