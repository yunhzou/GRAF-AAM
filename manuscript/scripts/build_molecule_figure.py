"""Vector chemical worked example; formal bond weights, no benchmark search.

RDKit generates chemical depictions. Its vector paths are transferred to
Matplotlib without rasterization. Analytical checks validate the illustrated
partial mappings, intrinsic response symmetry, and heavy-atom event lists.
"""
from pathlib import Path
import os,json,re,xml.etree.ElementTree as ET
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).resolve().parents[1]/"build/mpl-cache"))
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from matplotlib import pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch,Ellipse,FancyArrowPatch
from matplotlib.transforms import Affine2D

INK='#223340'; MUTED='#586974'; TEAL='#087F8C'; ORANGE='#CB702A'; PURPLE='#77559B'; RED='#B63F52'; BLUE='#3366A5'; LINE='#D5DEE2'
PALE_A=(.70,.88,.86); PALE_B=(1.,.82,.59); PALE_O=(.84,.76,.94)
ESTER='[CH3:1][C:2](=[O:3])[O:4][CH3:5]'
WATER='[OH2:6]'
ACID='[CH3:1][C:2](=[O:3])[OH:4]'
ALCOHOL='[CH3:5][OH:6]'
R=ESTER+'.'+WATER; P=ACID+'.'+ALCOHOL
ETHANOL='[CH3:11][CH2:12][OH:13]'
LACTATE='[CH3:21][CH:22]([OH:23])[C:24](=[O:25])[O:26][CH2:27][CH3:28]'

def weighted(smiles):
 m=Chem.MolFromSmiles(smiles); assert m is not None
 z={a.GetAtomMapNum():a.GetAtomicNum() for a in m.GetAtoms()}
 w={tuple(sorted((b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum()))):b.GetBondTypeAsDouble() for b in m.GetBonds()}
 return z,w

def check_example():
 zr,wr=weighted(R);zp,wp=weighted(P)
 identity={i:i for i in zr};exchange={**identity,3:4,4:3};challenger={**identity,4:6,6:4}
 def events(m):
  assert sorted(m.values())==sorted(zp) and all(zr[i]==zp[j] for i,j in m.items())
  return [[i,j,1 if d>0 else -1] for i in zr for j in zr if i<j and abs(d:=wp.get(tuple(sorted((m[i],m[j]))),0)-wr.get((i,j),0))>=.5]
 def preserves(m,fragments):
  return all(wp.get(tuple(sorted((m[i],m[j]))),0)>=.2 and abs(wp.get(tuple(sorted((m[i],m[j]))),0)-v)<=1 for f in fragments for (i,j),v in wr.items() if i in f and j in f)
 assert preserves(identity,[[1,2,3,4],[5],[6]])
 assert preserves(exchange,[[1,2,3,4],[5],[6]])
 assert preserves(challenger,[[1,2,3],[4,5],[6]])
 # The accepted acetoxy placement cannot grow through O4-C5 in either case.
 assert all(wp.get(tuple(sorted((m[4],m[5]))),0)<.2 for m in [identity,exchange])
 # The O3/O4 exchange is an exact symmetry of the matching-response graph,
 # although it is not a symmetry of the formal bond-order graph.
 response=lambda w:[w>=.2,*[w>=.2 and abs(w-v)<=1 for v in [1,2]]]
 assert response(1)==response(2)==[True,True,True]
 assert events(identity)==[[4,5,-1],[5,6,1]]
 assert events(exchange)==[[2,3,-1],[2,4,1],[4,5,-1],[5,6,1]]
 assert events(challenger)==[[2,4,-1],[2,6,1]]
 ze,we=weighted(ETHANOL);zl,wl=weighted(LACTATE)
 placements=[{11:21,12:22,13:23},{11:28,12:27,13:26}]
 for m in placements:
  assert len(set(m.values()))==3 and all(ze[i]==zl[j] for i,j in m.items())
  assert all(abs(wl.get(tuple(sorted((m[i],m[j]))),0)-v)<=1 and wl.get(tuple(sorted((m[i],m[j]))),0)>=.2 for (i,j),v in we.items())
 assert set(placements[0].values()) != set(placements[1].values())
 return {'placement_example':{'source_smiles':ETHANOL,'target_smiles':LACTATE,'two_admissible_images':placements,'scope':'Two selected placements, not an exhaustive placement list.'},'status':'passed','scope':'Constructed molecular illustration; heavy-atom projection with formal bond orders, not a measured trajectory or benchmark case. Hydrogens are hidden and are not counted here. Saved constraints still determine admissible full-family actions.',
 'reactants_smiles':R,'products_smiles':P,'iso_tol':1.,'edge_floor':.2,'event_threshold':.5,
 'mappings':{'retain_A':identity,'oxygen_exchange':exchange,'prioritize_B':challenger},
 'heavy_atom_events':{'retain_A':events(identity),'oxygen_exchange':events(exchange),'prioritize_B':events(challenger)},
 'response_labels':{'single':response(1),'double':response(2)},'saturated_fragment':[1,2,3,4]}

def svg_path(d):
 tokens=re.findall(r'[MLQCZmlqcz]|[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?',d)
 vertices=[];codes=[];i=0;start=(0,0)
 while i<len(tokens):
  c=tokens[i];i+=1
  assert c in {'M','L','Q','C','Z'},c
  if c=='Z':vertices.append(start);codes.append(MPath.CLOSEPOLY);continue
  n={'M':2,'L':2,'Q':4,'C':6}[c]
  a=list(map(float,tokens[i:i+n]));i+=n
  pts=list(zip(a[::2],a[1::2]));vertices.extend(pts)
  codes.extend([{'M':MPath.MOVETO,'L':MPath.LINETO,'Q':MPath.CURVE3,'C':MPath.CURVE4}[c]]*len(pts))
  if c=='M':start=pts[0]
 return MPath(vertices,codes)

def molecule(ax,smiles,x,y,w,h,owners=None,notes=None,font=24):
 m=Chem.MolFromSmiles(smiles);assert m is not None
 ids={a.GetAtomMapNum():a.GetIdx() for a in m.GetAtoms()}
 for a in m.GetAtoms():
  k=a.GetAtomMapNum();a.SetProp('atomNote',str((notes or {}).get(k,k)));a.SetAtomMapNum(0)
 rdDepictor.Compute2DCoords(m)
 # Shared acetyl coordinates make the ester and acid directly comparable.
 if {1,2,3,4}.issubset(ids):
  points={1:(-1.30,.75),2:(0,0),3:(0,-1.5),4:(1.30,.75),5:(2.60,0)}
  for k,i in ids.items():
   if k in points:m.GetConformer().SetAtomPosition(i,(*points[k],0))
 d=rdMolDraw2D.MolDraw2DSVG(int(w),int(h));o=d.drawOptions();o.clearBackground=False;o.useBWAtomPalette();o.fixedFontSize=font;o.annotationFontScale=.62;o.padding=.11;o.bondLineWidth=1.8
 colors={ids[k]:v for k,v in (owners or {}).items() if k in ids}
 bonds={b.GetIdx():colors[b.GetBeginAtomIdx()] for b in m.GetBonds() if b.GetBeginAtomIdx() in colors and colors.get(b.GetBeginAtomIdx())==colors.get(b.GetEndAtomIdx())}
 d.DrawMolecule(m,highlightAtoms=list(colors),highlightBonds=list(bonds),highlightAtomColors=colors,highlightBondColors=bonds)
 d.FinishDrawing();root=ET.fromstring(d.GetDrawingText())
 transform=Affine2D().translate(x,y)+ax.transData
 for e in root:
  kind=e.tag.split('}')[-1]
  if kind not in {'path','ellipse','rect'}:continue
  style=dict(item.split(':',1) for item in e.get('style','').split(';') if ':' in item)
  fill=style.get('fill',e.get('fill','none'));stroke=style.get('stroke',e.get('stroke','none'))
  kw=dict(facecolor=fill,edgecolor=stroke,linewidth=float(style.get('stroke-width','0').removesuffix('px'))*.60,transform=transform,zorder=3)
  if kind=='path':patch=PathPatch(svg_path(e.attrib['d']),**kw)
  elif kind=='ellipse':patch=Ellipse((float(e.attrib['cx']),float(e.attrib['cy'])),2*float(e.attrib['rx']),2*float(e.attrib['ry']),**kw)
  else:raise AssertionError('Unexpected background rectangle in molecule')
  ax.add_patch(patch)


def build(man):
 proof=check_example();(man/'evidence/molecule_example.json').write_text(json.dumps(proof,indent=2)+'\n')
 fig=plt.figure(figsize=(8.8,9.24));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1200),ylim=(1260,0));ax.axis('off')
 def t(x,y,s,size=9,color=INK,weight='normal',ha='left'):
  return ax.text(x,y,s,fontsize=size,color=color,weight=weight,ha=ha,va='center',linespacing=1.45)
 def ar(a,b,color=MUTED,rad=0,style='-|>'):
  ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=10,lw=1.1,color=color,connectionstyle=f'arc3,rad={rad}'))
 def head(x,y,letter,title,sub):
  t(x,y,letter,13,TEAL,'bold');t(x+34,y,title,11,INK,'bold');t(x,y+30,sub,8.6,MUTED)
 def rule(y):ax.plot([22,1178],[y,y],color=LINE,lw=.7)
 def acid(x,y,w,h,mapping=None,owners=None,letters=False):
  notes={i:chr(96+i) for i in range(1,7)} if letters else mapping
  molecule(ax,ACID,x,y,w,h,owners,notes)
 def products(x,y,scale=1.,mapping=None):
  acid(x,y,220*scale,136*scale,mapping)
  t(x+225*scale,y+65*scale,'+',11)
  molecule(ax,ALCOHOL,x+240*scale,y+24*scale,155*scale,85*scale,notes=mapping)

 head(24,24,'a','One reaction, several atom correspondences','Methyl acetate + water → acetic acid + methanol')
 molecule(ax,ESTER,20,68,345,136,font=25);t(377,130,'+',13)
 molecule(ax,WATER,405,90,110,82,font=25);ar((540,130),(650,130))
 acid(675,67,280,136,letters=True);t(970,130,'+',13)
 molecule(ax,ALCOHOL,994,91,185,89,notes={5:'e',6:'f'},font=25)
 t(265,218,'Reactant indices: 1–6',8.5,MUTED,ha='center');t(932,218,'Product indices: a–f',8.5,MUTED,ha='center')
 rule(242)
 head(24,273,'b','Growth and distinct placement branches','Local example: ethanol fragment in ethyl lactate.')
 for x,atoms,label in [(20,[11],'seed'),(203,[11,12],'extend'),(388,[11,12,13],'saturate')]:
  molecule(ax,ETHANOL,x,345,177,100,{i:PALE_A for i in atoms},{11:'u',12:'v',13:'w'},font=23)
  t(x+88,327,label,8.8,TEAL,'bold',ha='center')
 ar((181,393),(205,393));ar((368,393),(390,393))
 ar((466,439),(179,479),TEAL);ar((482,439),(443,479),TEAL)
 for x,atoms in [(15,[21,22,23]),(303,[26,27,28])]:
  molecule(ax,LACTATE,x,481,258,136,{i:PALE_A for i in atoms},{i:chr(97+i-21) for i in range(21,29)},font=23)
 t(148,627,'image {a,b,c}',9,ha='center');t(436,627,'image {f,g,h}',9,ha='center')
 t(24,655,'Different atom sets → different continuations (two shown).',8.4,MUTED)

 head(624,273,'c','Compress matching-preserving shuffles','Single and double C–O bonds both pass at τiso = 1.')
 acid(662,325,290,210,owners={3:PALE_O,4:PALE_O},letters=True)
 t(986,364,'R bond     P bond',8.6,INK,'bold')
 t(986,401,'2  →  2 or 1',9);t(986,436,'1  →  1 or 2',9)
 t(986,481,'|Δw| ≤ 1',10,TEAL,'bold')
 t(644,558,'g = (c d)',12,PURPLE,'bold');t(819,558,'one coupled O exchange',9,PURPLE)
 t(644,597,'χ(1) = χ(2): an intrinsic fragment symmetry.',9)
 t(644,638,'Saved constraints decide which shuffles are allowed.\nThis is not formal-bond-order chemical symmetry.',8.5,MUTED)
 rule(683)

 head(24,715,'d','Branch again when fragments compete','Retain A, or let B claim the contacting oxygen 4.')
 molecule(ax,ESTER,25,763,310,153,{**{i:PALE_A for i in [1,2,3,4]},5:PALE_B},font=25)
 t(365,802,'A = {1,2,3,4}',9,TEAL,'bold');t(365,837,'B starts at {5}',9,ORANGE,'bold')
 ar((248,881),(111,927),TEAL);ar((310,881),(411,927),ORANGE)
 t(128,931,'keep A',9.5,TEAL,'bold',ha='center');t(417,931,'prioritize B',9.5,ORANGE,'bold',ha='center')
 molecule(ax,ESTER,16,949,253,133,{**{i:PALE_A for i in [1,2,3,4]},5:PALE_B},font=23)
 molecule(ax,ESTER,306,949,253,133,{**{i:PALE_A for i in [1,2,3]},4:PALE_B,5:PALE_B},font=23)
 t(138,1096,'4 → d;  5 → e',9,ha='center');t(431,1096,'4 → f;  5 → e',9,ha='center')
 t(138,1124,'cleave O4–C5',8.5,TEAL,ha='center');t(431,1124,'cleave C2–O4',8.5,ORANGE,ha='center')

 head(624,715,'e','Decode different bond-event patterns','Products labeled by their assigned reactant indices.')
 products(640,754,.99)
 t(1086,807,'2',17,TEAL,'bold',ha='center');t(1086,839,'events',8.5,TEAL,ha='center')
 t(645,908,'−(4,5), +(5,6)',9,TEAL)
 products(640,924,.99,mapping={1:1,2:2,3:4,4:3,5:5,6:6})
 t(1086,977,'4',17,RED,'bold',ha='center');t(1086,1009,'events',8.5,RED,ha='center')
 t(645,1078,'Same changes, plus −(2,3), +(2,4)',9,RED)
 t(645,1120,'O exchange changes events despite admissible matching.',8.3,MUTED)
 rule(1144)
 head(24,1174,'f','One unordered branch; retain its mapping families','The oxygen assignments in (c,e) share {(1,2,3,4) ↔ (a,b,c,d), 5 ↔ e, 6 ↔ f}, regardless of growth order.')
 t(24,1240,'Constructed heavy-atom illustration · formal bond weights · event threshold 0.5 · H transfers omitted',8.4,MUTED)
 for ext in ['pdf','svg','png']:fig.savefig(man/'figs'/f'fig1_algorithm.{ext}',facecolor='white',dpi=260)
 plt.close(fig)
 print('Molecular figure: RDKit vector structures; all example assertions passed.')

if __name__=='__main__':build(Path(__file__).resolve().parents[1])
