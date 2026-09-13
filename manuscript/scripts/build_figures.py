"""Build the final paper's vector figures from checked-in evidence; no AAM runs."""
import os,json
from pathlib import Path
MAN=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(MAN/'build/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch,FancyBboxPatch
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':240})
INK='#223340';MUTED='#586974';GREEN='#087F8C';BLUE='#3366A5';ORANGE='#CB702A';RED='#B63F52';LIGHT='#D5DEE2';PURPLE='#77559B'
E=MAN/'evidence';F=MAN/'figs';F.mkdir(exist_ok=True)
def read(n):return json.loads((E/n).read_text())
def save(fig,n):
 for ext in ['pdf','svg','png']:fig.savefig(F/f'{n}.{ext}',facecolor='white',bbox_inches='tight',pad_inches=.05)
 plt.close(fig)
def text(ax,x,y,s,**kw):
 opts=dict(color=INK,fontsize=9,ha='left',va='center');opts.update(kw);opts['fontsize']=max(9,opts['fontsize']);return ax.text(x,y,s,**opts)
def arrow(ax,a,b,color=MUTED,style='-|>',rad=0,lw=1.3):
 ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=10,color=color,lw=lw,connectionstyle=f'arc3,rad={rad}'))
def node(ax,x,y,label='',color=INK,fill='white',r=.14):
 ax.add_patch(Circle((x,y),r,facecolor=fill,edgecolor=color,lw=1.4,zorder=4))
 text(ax,x,y,label,color=color,ha='center',fontsize=8,zorder=5)
def edge(ax,a,b,color=LIGHT,lw=2,ls='-'):
 ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=lw,ls=ls,zorder=1)
def panel(ax,letter,title,subtitle):
 ax.set(xlim=(0,5),ylim=(0,3.15));ax.axis('off')
 text(ax,.05,3.0,letter,weight='bold',fontsize=13,color=GREEN)
 text(ax,.34,3.0,title,weight='bold',fontsize=11)
 text(ax,.05,2.70,subtitle,fontsize=8.4,color=MUTED)
 ax.plot([.05,4.92],[2.51,2.51],color=LIGHT,lw=.7)
def tag(ax,x,y,w,s,color=GREEN):
 ax.add_patch(FancyBboxPatch((x,y-.17),w,.34,boxstyle='round,pad=.04,rounding_size=.05',fc='white',ec=color,lw=1))
 text(ax,x+w/2,y,s,ha='center',color=color,fontsize=8)
def pairgraph(ax,x,y,labels=('C','O'),weight=None,color=GREEN):
 edge(ax,(x,y),(x+.85,y),color)
 node(ax,x,y,labels[0],color);node(ax,x+.85,y,labels[1],color)
 if weight is not None:text(ax,x+.425,y+.28,weight,ha='center',fontsize=8)
def vgraph(ax,x,y,weights,labels=('C','O$_1$','O$_2$'),color=GREEN):
 p=[(x,y),(x-.63,y+.7),(x+.63,y+.7)]
 for i in [1,2]:edge(ax,p[0],p[i],color,1.6)
 for xy,label in zip(p,labels):node(ax,*xy,label,color=color,r=.16)
 for i,w in zip([1,2],weights):
  text(ax,x+(-.52 if i==1 else .52),y+.24,w,ha='center',fontsize=8)

fig,axs=plt.subplots(3,2,figsize=(8.4,9.0))
fig.subplots_adjust(left=.025,right=.985,bottom=.025,top=.985,wspace=.10,hspace=.16)
a=axs[0,0];panel(a,'a','Weighted matching','R–P weight compatibility controls admissible placements.')
text(a,.33,2.24,'Reactant edge',fontsize=8.5,weight='bold')
pairgraph(a,.48,1.81,weight='1.20')
arrow(a,(1.72,1.81),(2.27,2.00));arrow(a,(1.72,1.76),(2.27,1.00))
pairgraph(a,2.55,2.00,weight='0.90');text(a,3.70,2.00,'accept',color=GREEN,fontsize=8.5)
pairgraph(a,2.55,1.00,weight='2.40',color=RED);text(a,3.70,1.00,'reject',color=RED,fontsize=8.5)
text(a,2.65,1.54,r'$|1.20-0.90|=0.30$',fontsize=8.5)
text(a,2.65,.56,r'$|1.20-2.40|=1.20$',fontsize=8.5)
text(a,.08,.22,r'Example: $\tau_{\mathrm{iso}}=1.0$; edge floor $b=0.2$.',fontsize=8.3)

b=axs[0,1];panel(b,'b','Grow, refine, then branch','Grow even when only one placement remains.')
for k,center in enumerate([.65,2.25,3.85]):
 xs=[center-.39,center,center+.39]
 for i in range(2):edge(b,(xs[i],1.86),(xs[i+1],1.86),GREEN if i<k else LIGHT)
 for i,x in enumerate(xs):node(b,x,1.86,'',GREEN if i<=k else LIGHT,fill=GREEN if i<=k else 'white',r=.09)
 text(b,center,2.23,['seed','extend','saturate'][k],ha='center',fontsize=8.5,weight='bold')
 text(b,center,1.50,['placements','refined set','survivors'][k],ha='center',fontsize=7.8)
 for j in range([3,2,2][k]):tag(b,center-.42,1.18-j*.40,.84,f'P{j+1}',BLUE)
 if k<2:arrow(b,(center+.55,1.86),(center+1.01,1.86))
arrow(b,(4.33,1.18),(4.79,1.34),BLUE);arrow(b,(4.33,.78),(4.79,.54),BLUE)
text(b,.08,-.04,'Continue remaining atoms under each placement.',fontsize=8.3)

c=axs[1,0];panel(c,'c','Compress with automorphisms','Generators preserve matching tests, not numerical weights.')
vgraph(c,1.04,1.14,['1.00','0.90']);vgraph(c,3.68,1.14,['0.51','0.49'],labels=('C','O$_a$','O$_b$'),color=BLUE)
text(c,1.04,2.29,'R fragment',ha='center',fontsize=8.5,weight='bold');text(c,3.68,2.29,'P fragment',ha='center',fontsize=8.5,weight='bold')
arrow(c,(1.89,1.57),(2.78,1.57));text(c,2.34,1.97,r'$\tau_{\mathrm{iso}}=1$',ha='center',fontsize=8.2)
arrow(c,(3.10,2.00),(4.24,2.00),PURPLE,style='<->',rad=-.20)
text(c,.09,.64,r'Generator: $g=(O_a\ O_b)$',color=PURPLE,fontsize=9)
text(c,.09,.27,'Both O assignments pass the matching tests.\nThe group encodes coupled atom moves.',fontsize=8.3,linespacing=1.55)

d=axs[1,1];panel(d,'d','Let neighboring fragments compete','Keep the original branch; add a challenger continuation.')
def chain(y,owners,label):
 xs=np.linspace(.73,2.66,5)
 for i in range(4):edge(d,(xs[i],y),(xs[i+1],y),GREEN if owners[i]==owners[i+1]=='A' else ORANGE if owners[i]==owners[i+1]=='B' else MUTED,lw=2,ls='-' if owners[i]==owners[i+1] else ':')
 for i,(x,owner) in enumerate(zip(xs,owners)):node(d,x,y,str(i+1),GREEN if owner=='A' else ORANGE,r=.12)
 text(d,3.06,y,label,fontsize=8.3)
chain(2.13,'AAABB','A holds its assignments')
arrow(d,(1.69,1.88),(1.06,1.39),GREEN);arrow(d,(2.18,1.88),(2.68,1.39),ORANGE)
text(d,.49,1.29,'keep A',color=GREEN,fontsize=8.4,weight='bold');text(d,2.09,1.29,'prioritize B',color=ORANGE,fontsize=8.4,weight='bold')
# Two branches as miniature partition strips.
for offset,owners in [(.17,'AAABB'),(2.76,'AABBB')]:
 for i,owner in enumerate(owners):node(d,offset+.18+i*.37,.88,str(i+1),GREEN if owner=='A' else ORANGE,r=.11)
text(d,.05,.24,'Rematch B at the contact; release conflicts,\nanchor unaffected atoms, and regrow gaps.',fontsize=8.3,linespacing=1.55)

e=axs[2,0];panel(e,'e','Deduplicate final fragment pairs','Growth order disappears; internal alternatives remain.')
for y,names in [(2.12,['A','B']),(1.49,['B','A'])]:
 for k,name in enumerate(names):tag(e,.12+k*.96,y,.61,name,GREEN if name=='A' else ORANGE)
 arrow(e,(.79,y),(.99,y))
arrow(e,(1.77,2.12),(2.29,1.87));arrow(e,(1.77,1.49),(2.29,1.76))
tag(e,2.36,1.84,2.34,r'$\{(R_A,P_A),(R_B,P_B)\}$',BLUE)
text(e,3.53,1.37,'one final branch',ha='center',color=BLUE,fontsize=9,weight='bold')
text(e,.10,.80,r'Branch contents: $\mathcal{F}_1\,\cup\,\mathcal{F}_2\,\cup\,\cdots$',fontsize=10)
text(e,.10,.34,'Share identical records; retain different\nconstraints, actions, and witnesses.',fontsize=8.3,linespacing=1.55)

f=axs[2,1];panel(f,'f','Decode distinct bond-event classes','Project allowed mappings onto signed bond changes.')
text(f,.09,2.20,'O assignment in (c)',fontsize=8.6,weight='bold');text(f,2.63,2.20,r'$\Delta W$',fontsize=8.6,weight='bold');text(f,4.08,2.20,'events',fontsize=8.6,weight='bold')
text(f,.09,1.81,r'$O_1\!\mapsto O_a,\quad O_2\!\mapsto O_b$',fontsize=9)
text(f,2.60,1.81,'−0.49, −0.41',fontsize=8.4);text(f,4.34,1.81,'0',color=GREEN,weight='bold',fontsize=11)
text(f,.09,1.37,r'$O_1\!\mapsto O_b,\quad O_2\!\mapsto O_a$',fontsize=9)
text(f,2.60,1.37,'−0.51, −0.39',fontsize=8.4);text(f,4.34,1.37,'1',color=RED,weight='bold',fontsize=11)
text(f,.09,.94,r'At $\delta=0.5$: the second mapping weakens $C{-}O_1$.',fontsize=8.2)
for x,w,s in [(.12,1.38,'find new class'),(1.86,1.18,'save witness'),(3.41,1.28,'exclude class')]:tag(f,x,.51,w,s,BLUE)
arrow(f,(1.57,.51),(1.78,.51),BLUE);arrow(f,(3.10,.51),(3.34,.51),BLUE)
text(f,.10,.09,'Repeat until the event window is exhausted.',fontsize=8.2)
# The schematic's two mappings must produce the claimed event counts.
r=np.array([1.,.9]);p=np.array([.51,.49]);assert [int((np.abs(q-r)>=.5).sum()) for q in [p,p[::-1]]]==[0,1]
assert np.all(np.abs(p[:,None]-r[None,:])<=1.)
save(fig,'fig1_algorithm')

seed=read('seed_comparison.json');methods=seed['methods'];slap=read('slap_sweep.json');N=1851
keys=['seeds1','seeds2','seeds3','seeds10'];counts=[methods[k]['golden_outcomes']['recovered'] for k in keys];costs=[methods[k]['common_mean_cpu_seconds'] for k in keys]
assert seed['fresh'] and all(sum(methods[k]['golden_outcomes'].values()) == N for k in keys)
assert len(seed['common_case_indices']) > 0
fig,axs=plt.subplots(1,2,figsize=(8.4,3.55),gridspec_kw={'width_ratios':[1.45,1]},layout='constrained')
a=axs[0];labels=['GRAFT · 1 order','GRAFT · 10 orders','SLAP sweep'];vals=[counts[0],counts[3],slap['sweep_union_recovered']];y=np.arange(3)
a.barh(y,np.array(vals)/N*100,color=[GREEN,BLUE,ORANGE],height=.50)
for i,v in enumerate(vals):a.text(v/N*100-1.5,i,f'{v:,} / {N:,}  ({v/N*100:.2f}%)',ha='right',va='center',color='white',fontsize=9,weight='bold')
a.set(yticks=y,yticklabels=labels,xlim=(0,102),xticks=[0,25,50,75,100],xlabel='Verified reference-family recovery (%)');a.invert_yaxis();a.set_title('a  Golden mapping coverage',loc='left',weight='bold',pad=13)
a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
a=axs[1];shown=[costs[0],costs[1],seed['same_host_timing']['metrics']['slap']['mean_seconds']];y=np.arange(3);a.barh(y,shown,color=[GREEN,BLUE,ORANGE],height=.53)
for i,v in enumerate(shown):a.text(v+.10,i,f'{v:.2f}',va='center',fontsize=9)
a.set(yticks=y,yticklabels=['GRAFT · 1 order','GRAFT · 2 orders','SLAP sweep'],xlim=(0,max(shown)*1.23),xlabel='Mean CPU seconds / reaction');a.invert_yaxis();a.set_title('b  Same-Mac workflow cost',loc='left',weight='bold',pad=13);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
save(fig,'fig2_golden')
table=[]
for k,label in zip(keys,['1 seed order','2 seed orders','3 seed orders','10 seed orders']):
 d=methods[k];o=d['golden_outcomes'];cost_text=f"{d['common_mean_cpu_seconds']:.2f}" if d['common_mean_cpu_seconds'] is not None else '--';table.append(f"{label} & {o['recovered']:,} & {d['golden_recovery_percent']:.2f} & {o['not_recovered']} & {o['unknown']} & {cost_text} \\\\")
(MAN/'includes/generated-seed-table.tex').write_text('\\begin{tabular}{lrrrrr}\\toprule\nConfiguration & Recovered & \\% & Absent & Unknown & CPU s/reaction\\\\\\midrule\n'+'\n'.join(table)+'\n\\bottomrule\\end{tabular}\n')

comp=read('competition_final.json');dedup=read('final_dedup.json');rows=dedup['per_case']
assert len(rows)==140 and sum(r['final_branches'] for r in rows)==dedup['new_branches']
assert dedup['fresh_decode_all140'] and all(r['complete'] for r in rows)
fig,axs=plt.subplots(1,2,figsize=(8.4,3.8),gridspec_kw={'width_ratios':[1.1,1.2]},layout='constrained')
a=axs[0];names=['SLAP sweep','Native SLAP'];cs=[comp['comparisons'][k] for k in ['slap_sweep','native_slap']]
for i,d in enumerate(cs):
 n=d['union_classes'];total=d['total_classes'];a.barh(i,n,color=GREEN,height=.5);a.barh(i,total-n,left=n,color=RED,height=.5)
 a.text(n/2,i,f'{n}/{total} classes recovered',ha='center',va='center',color='white',fontsize=9,weight='bold')
a.set(yticks=[0,1],yticklabels=names,xlim=(0,180),xticks=[0,50,100,150],xlabel='Minimum-event classes in comparator output');a.invert_yaxis();a.set_title('a  Coverage of comparator alternatives',loc='left',weight='bold',pad=13);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
a.text(.02,.05,'All comparator classes covered in 139/140 reactions.\nRed: missing classes, all in one reaction.',transform=a.transAxes,fontsize=8.4,color=MUTED,va='bottom');a.set_ylim(1.85,-.6)
a=axs[1];old=np.array([r['old_literal_branches'] for r in rows]);new=np.array([r['final_branches'] for r in rows]);a.plot([10,8e4],[10,8e4],color=LIGHT,ls='--',lw=1,label='equal counts');a.scatter(old,new,s=21,c=GREEN,alpha=.65,edgecolors='white',linewidths=.4)
r=next(r for r in rows if r['case']==25);a.scatter([r['old_literal_branches']],[r['final_branches']],s=42,c=ORANGE,zorder=5);a.annotate('52,669 → 2,856',xy=(52669,2856),xytext=(1000,130),arrowprops=dict(arrowstyle='-',color=ORANGE),color=ORANGE,fontsize=8.5)
a.set(xscale='log',yscale='log',xlim=(12,8e4),ylim=(12,8e4),xlabel='Ordered branches per reaction',ylabel='Unordered final branches per reaction');a.set_title('b  Final fragment deduplication',loc='left',weight='bold',pad=13);a.grid(alpha=.12)
a.text(.03,.97,'Total: 237,645 → 124,641\nMedian: 596.5 → 397',transform=a.transAxes,va='top',fontsize=8.5)
save(fig,'fig3_coordinate')
windows=sorted((int(k),v) for k,v in dedup['window_distribution'].items())
(MAN/'includes/generated-decoder-table.tex').write_text('\\begin{tabular}{l'+ 'r'*len(windows)+'}\\toprule\nMaximum events & '+' & '.join(str(k) for k,v in windows)+r'\\'+'\nReactions (all complete) & '+' & '.join(str(v) for k,v in windows)+r'\\'+'\n'+r'\bottomrule\end{tabular}'+'\n')
paired=slap['paired_recovery']['seeds1']
shuffle=[len(r['new_classes_absent_from_repair_representatives']) for r in comp['cases']]
assert sum(v>0 for v in shuffle)==1
assert sum(shuffle)==1
macros=dict(GoldenOneRecovered=f"{counts[0]:,}",GoldenOnePercent=f"{100*counts[0]/N:.2f}",
 GoldenThreeRecovered=f"{counts[2]:,}",GoldenThreePercent=f"{100*counts[2]/N:.2f}",
 GoldenTenRecovered=f"{counts[3]:,}",GoldenTenPercent=f"{100*counts[3]/N:.2f}",GoldenTenUnknown=str(methods["seeds10"]["golden_outcomes"]["unknown"]),SlapRecovered=f"{slap['sweep_union_recovered']:,}",
 SlapPercent=f"{100*slap['sweep_union_recovered']/N:.2f}",
 RecoveryDifference=f"{100*(counts[0]-slap['sweep_union_recovered'])/N:.2f}",
 CommonCases=f"{len(seed['common_case_indices']):,}",OneCPU=f"{costs[0]:.2f}",TwoCPU=f"{costs[1]:.2f}",SlapCPU=f"{shown[2]:.2f}",
 AamOnly=str(len(paired['aam_only'])),SlapOnly=str(len(paired['slap_only'])),
 SweepRecovered=str(comp['comparisons']['slap_sweep']['union_classes']),
 SweepClasses=str(comp['comparisons']['slap_sweep']['total_classes']),
 CompetitionClasses=str(comp['new_window_class_count']),CompetitionCases=str(len(comp['new_window_cases'])),
 ShuffleClasses=str(sum(shuffle)),DecodeWallMinutes=f"{dedup['wall_seconds']/60:.2f}",
 DecodeCPUMinutes=f"{dedup['cpu_seconds']/60:.2f}",DecodeGiB=f"{dedup['peak_mib']/1024:.2f}",
 SlapAbsent=str(slap['outcomes']['not_recovered']),SlapUnknown=str(slap['outcomes']['unknown']),
 FlatFamilies=f"{dedup['flat_families']:,}",DecodedClasses=str(dedup['event_classes']))
(MAN/'includes/generated-results.tex').write_text('% Generated from complete fresh benchmark evidence.\n'+''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in macros.items()))
print('Built three main figures and complete-campaign tables and numerical macros.')
