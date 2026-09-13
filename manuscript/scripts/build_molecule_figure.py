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
from matplotlib.patches import PathPatch,Ellipse,FancyArrowPatch,FancyBboxPatch,Circle
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
  k=a.GetAtomMapNum()
  if notes is None:a.SetProp('atomNote',str(k))
  elif notes is not False and k in notes:a.SetProp('atomNote',str(notes[k]))
  a.SetAtomMapNum(0)
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
 return {k:(x+d.GetDrawCoords(i).x,y+d.GetDrawCoords(i).y) for k,i in ids.items()}


def build(man):
 proof=check_example();(man/'evidence/molecule_example.json').write_text(json.dumps(proof,indent=2)+'\n')
 # A connected branching overview, with chemical drawings as the tree states.
 fig=plt.figure(figsize=(10.8,7.25));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1440),ylim=(967,0));ax.axis('off')
 navy='#22324D';muted='#64748B';faint='#BCC7D5';bg='#F5F7FB';violet='#7957BA';orange='#CC762F';green='#128477'
 def t(x,y,s,size=11,color=navy,weight='normal',ha='left'):
  return ax.text(x,y,s,fontsize=size,color=color,weight=weight,ha=ha,va='center',linespacing=1.3,zorder=6)
 def box(x,y,w,h,fc='white',ec=LINE,lw=.8,r=12):
  patch=FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',facecolor=fc,edgecolor=ec,lw=lw,zorder=1)
  ax.add_patch(patch)
 def ar(a,b,color=navy,lw=1.9,dashed=False,rad=0,style='-|>'):
  ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=11,lw=lw,color=color,linestyle=(0,(3,3)) if dashed else '-',connectionstyle=f'arc3,rad={rad}',zorder=2))
 def tree(points,color=navy,lw=2.0):
  xs,ys=zip(*points[:-1]);ax.plot(xs,ys,color=color,lw=lw,solid_joinstyle='round',solid_capstyle='round',zorder=2)
  ar(points[-2],points[-1],color,lw)
 def dot(x,y,color=navy,r=4.3):ax.add_patch(Circle((x,y),r,fc=color,ec='none',zorder=3))
 def label(x,y,letter,title):
  t(x,y,letter,14,navy,'bold');t(x+36,y,title,14,navy,'bold')
 def continuation(x,y):
  tree([(x,y),(x+35,y),(x+35,y-27),(x+96,y-27)],faint,1.3)
  tree([(x,y),(x+35,y),(x+35,y+27),(x+96,y+27)],faint,1.3)
  dot(x+35,y,faint,3)
  for yy in [y-27,y+27]:t(x+114,yy,'···',13,faint)

 # A: the same growing ethanol fragment has distinct target atom-set images.
 box(12,12,1416,437,bg,bg,r=16)
 label(34,47,'a','Grow a fragment, then branch over its surviving placements')
 t(34,79,'Local matching example: ethanol in ethyl lactate',10.5,muted)
 for x,atoms,name in [(34,[11],'Seed'),(220,[11,12],'Extend'),(406,[11,12,13],'Saturate')]:
  t(x+75,140,name,11,green,'bold',ha='center')
  molecule(ax,ETHANOL,x,175,155,132,{i:PALE_A for i in atoms},notes=False,font=32)
 ar((185,242),(214,242),green,1.8);ar((371,242),(400,242),green,1.8)
 t(299,344,'Grow while refining the placement set',10.5,muted,ha='center')
 t(299,383,r'$|w_R-w_P|\leq\tau_{\mathrm{iso}}$',12,navy,ha='center')
 ar((565,242),(630,242),navy,2.2)
 tree([(630,242),(667,242),(667,170),(731,170)],navy,2.2)
 tree([(630,242),(667,242),(667,357),(731,357)],navy,2.2)
 dot(667,242)
 for x,y,atoms,name in [(738,105,[21,22,23],'Placement 1'),(738,292,[26,27,28],'Placement 2')]:
  box(x,y,345,142,'white','#D6E1E9',.9,10)
  molecule(ax,LACTATE,x+15,y+4,244,134,{i:PALE_A for i in atoms},notes=False,font=31)
  t(x+256,y+53,name.replace(' ','\n'),10.5,navy,'bold',ha='center')
  continuation(x+345,y+71)
 t(1265,253,'Continue\nremaining atoms',10.5,muted,ha='center')
 t(1144,79,'Different target atom sets',10.5,navy,'bold',ha='center')

 # B: the existing continuation survives; a challenger makes a new branch.
 box(12,467,702,448,bg,bg,r=16)
 label(34,503,'b','Add a branch when fragments compete')
 t(34,536,'Methyl acetate hydrolysis',10.5,muted)
 molecule(ax,ESTER,22,642,260,162,{**{i:PALE_A for i in [1,2,3,4]},5:PALE_B},notes={4:'4'},font=31)
 t(149,819,r'$+\ \mathrm{H_2O}$',13,navy,ha='center')
 t(151,596,'A',12,green,'bold',ha='center');t(192,596,'B',12,orange,'bold',ha='center')
 t(150,859,'Contact at O4',10.5,muted,ha='center')
 ar((280,722),(301,722),navy,2.2)
 tree([(301,722),(326,722),(326,628),(379,628)],green,2.2)
 tree([(301,722),(326,722),(326,818),(379,818)],orange,2.2)
 dot(326,722)
 for y,owners,title,color in [(560,{**{i:PALE_A for i in [1,2,3,4]},5:PALE_B},'Keep A',green),(750,{**{i:PALE_A for i in [1,2,3]},4:PALE_B,5:PALE_B},'Give B priority',orange)]:
  box(385,y,264,136,'white','#D6E1E9',.9,10)
  t(518,y-16,title,11,color,'bold',ha='center')
  molecule(ax,ESTER,391,y+2,250,130,owners,notes={4:'4'},font=30)
  if title=='Keep A':tree([(650,y+68),(729,y+68),(729,714),(753,714)],color,1.8)
  else:
   ar((650,y+68),(685,y+68),color,1.8)
   t(695,y+68,'···',12,color)
 # The orange path is an added continuation, not replacement of the green path.
 t(370,901,'Keep the original path; add the challenger.',10.5,muted,ha='center')

 # C: output projection uses dashed arrows; it is not a search branch.
 box(732,467,696,448,bg,bg,r=16)
 label(754,503,'c','Decode each compressed family')
 t(754,536,'The oxygen exchange preserves matching at '+r'$\tau_{\mathrm{iso}}=1$',10.5,muted)
 box(757,621,216,189,'white','#CEC1E3',1.0,12)
 coords=molecule(ax,ACID,759,630,210,159,{3:PALE_O,4:PALE_O},notes={3:'c',4:'d'},font=33)
 t(865,594,'Matching symmetry',10.5,violet,'bold',ha='center')
 t(865,790,r'$g=(c\ d)$',12,violet,ha='center')
 t(865,851,'One branch\ncorrelated alternatives',10.5,violet,ha='center')
 # Witness cards show oxygen identities while leaving unchanged carbons unlabeled.
 for y,mapping,count,color in [(575,{3:3,4:4,5:5,6:6},2,green),(757,{3:4,4:3,5:5,6:6},4,RED)]:
  box(1056,y,347,144,'white','#D6E1E9',.9,10)
  molecule(ax,ACID,1060,y+10,146,118,notes={i:mapping[i] for i in [3,4]},font=29)
  t(1209,y+60,'+',11)
  molecule(ax,ALCOHOL,1220,y+39,95,72,notes={6:'6'},font=28)
  t(1359,y+51,str(count),20,color,'bold',ha='center')
  t(1359,y+83,'events',10,color,ha='center')
 ar((975,704),(1053,647),violet,1.8,True)
 ar((975,730),(1053,829),violet,1.8,True)
 t(1230,739,'Exchange O3 / O4',10,violet,ha='center')
 # Shared legend, kept outside the chemical trees.
 ar((35,949),(84,949),navy,1.7);t(95,949,'Search continuation',10,muted)
 ar((394,949),(443,949),violet,1.7,True);t(454,949,'Event projection',10,muted)
 t(1407,949,'Formal bond orders · heavy-atom events at δ = 0.5',10,muted,ha='right')
 fig.canvas.draw()
 renderer=fig.canvas.get_renderer()
 for label in ax.texts:
  extent=label.get_window_extent(renderer)
  assert extent.x0>=-1 and extent.x1<=fig.bbox.width+1 and extent.y0>=-1 and extent.y1<=fig.bbox.height+1,label.get_text()
 for ext in ['pdf','svg','png']:fig.savefig(man/'figs'/f'fig1_algorithm.{ext}',facecolor='white',dpi=260)
 plt.close(fig)
 print('Branching overview: RDKit vector molecules and checked illustrative mappings.')

if __name__=='__main__':build(Path(__file__).resolve().parents[1])
