"""Build the final paper's vector figures from checked-in evidence; no AAM runs."""
import os,json
from pathlib import Path
MAN=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(MAN/'build/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':240})
INK='#223340';MUTED='#586974';GREEN='#087F8C';BLUE='#3366A5';ORANGE='#CB702A';RED='#B63F52';LIGHT='#D5DEE2';PURPLE='#77559B'
E=MAN/'evidence';F=MAN/'figs';F.mkdir(exist_ok=True)
def read(n):return json.loads((E/n).read_text())
def save(fig,n):
 for ext in ['pdf','svg','png']:fig.savefig(F/f'{n}.{ext}',facecolor='white',bbox_inches='tight',pad_inches=.05)
 plt.close(fig)
from build_molecule_figure import build as build_molecule_figure
build_molecule_figure(MAN)

seed=read('seed_comparison.json');methods=seed['methods'];slap=read('slap_sweep.json');N=1851
keys=['seeds1','seeds2','seeds3','seeds10'];counts=[methods[k]['golden_outcomes']['recovered'] for k in keys];costs=[methods[k]['common_mean_cpu_seconds'] for k in keys]
assert seed['fresh'] and all(sum(methods[k]['golden_outcomes'].values()) == N for k in keys)
assert len(seed['common_case_indices']) > 0
fig,axs=plt.subplots(1,2,figsize=(8.4,3.8),gridspec_kw={'width_ratios':[1.45,1]},layout='constrained')
a=axs[0];labels=['GRAFT · 1 seed\n(default)','GRAFT · 2 seeds','GRAFT · 3 seeds','GRAFT · 10 seeds','SLAP sweep'];vals=[*counts,slap['sweep_union_recovered']];y=np.arange(5)
a.barh(y,np.array(vals)/N*100,color=[GREEN,BLUE,'#52789F',PURPLE,ORANGE],height=.58)
for i,v in enumerate(vals):a.text(v/N*100-1.5,i,f'{v:,} / {N:,}  ({v/N*100:.2f}%)',ha='right',va='center',color='white',fontsize=9,weight='bold')
a.set(yticks=y,yticklabels=labels,xlim=(0,102),xticks=[0,25,50,75,100],xlabel='Verified reference-family recovery (%)');a.invert_yaxis();a.set_title('a  Recovery with sweeps',loc='left',weight='bold',pad=13)
a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
a=axs[1];shown=[costs[0],costs[1],seed['same_host_timing']['metrics']['slap']['mean_seconds']];y=np.arange(3);a.barh(y,shown,color=[GREEN,BLUE,ORANGE],height=.53)
for i,v in enumerate(shown):a.text(v+.10,i,f'{v:.2f}',va='center',fontsize=9)
a.set(yticks=y,yticklabels=['GRAFT · 1 seed\n(default)','GRAFT · 2 seeds','SLAP sweep'],xlim=(0,max(shown)*1.23),xlabel='Mean CPU seconds / reaction');a.invert_yaxis();a.set_title('b  Sweep workflow cost (Mac)',loc='left',weight='bold',pad=13);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
save(fig,'fig2_golden')
uncut=read('unswept.json');table=[]
for name,label in [('graft','GRAFT, no sweep (ablation)'),('slap','SLAP, no sweep')]:
 o=uncut['methods'][name]['counts']
 table.append(f"{label} & {o['recovered']:,} & {100*o['recovered']/N:.2f} & {o['not_recovered']} & {o.get('unknown',0)} & -- " + r"\\")
ablations=table;table=[]
for k,label in zip(keys,['GRAFT, 1 seed + sweep (default)','GRAFT, 2 seeds + sweep','GRAFT, 3 seeds + sweep','GRAFT, 10 seeds + sweep']):
 d=methods[k];o=d['golden_outcomes'];cost_text=f"{d['common_mean_cpu_seconds']:.2f}" if d['common_mean_cpu_seconds'] is not None else '--';table.append(f"{label} & {o['recovered']:,} & {d['golden_recovery_percent']:.2f} & {o['not_recovered']} & {o['unknown']} & {cost_text} \\\\")
table.append(f"SLAP sweep & {slap['outcomes']['recovered']:,} & {100*slap['outcomes']['recovered']/N:.2f} & {slap['outcomes']['not_recovered']} & {slap['outcomes']['unknown']} & {shown[2]:.2f} " + r"\\")
table += [r'\midrule'] + ablations
(MAN/'includes/generated-seed-table.tex').write_text('\\begin{tabular}{lrrrrr}\\toprule\nConfiguration & Recovered & \\% & Absent & Unknown & CPU s/reaction\\\\\\midrule\n'+'\n'.join(table)+'\n\\bottomrule\\end{tabular}\n')

comp=read('competition_final.json');dedup=read('final_dedup.json');rows=dedup['per_case']
assert len(rows)==140 and sum(r['final_branches'] for r in rows)==dedup['new_branches']
assert dedup['fresh_decode_all140'] and all(r['complete'] for r in rows)
fig,axs=plt.subplots(2,2,figsize=(8.4,6.0),layout='constrained',gridspec_kw={'hspace':.17,'wspace':.10})
a=axs[0,0];total=dedup['event_classes'];added=comp['new_window_class_count'];base=total-added
assert added==sum(len(r['new_window_classes']) for r in comp['cases'])
a.barh(0,base,color=GREEN,height=.48,label='Ordinary search')
a.barh(0,added,left=base,color=PURPLE,height=.48,label='Added by competition')
a.text(base/2,0,str(base),ha='center',va='center',color='white',weight='bold')
a.text(base+added/2,0,str(added),ha='center',va='center',color='white',weight='bold')
a.text(.02,.94,f'{total} event classes · all 140 windows complete',transform=a.transAxes,va='top',weight='bold',fontsize=9)
a.text(.02,.12,f'Competition adds {added} classes in {len(comp["new_window_cases"])} reactions.',transform=a.transAxes,fontsize=8.4,color=MUTED)
a.set(yticks=[],xlim=(0,360),ylim=(-.85,.8),xticks=[0,100,200,300],xlabel='GRAFT classes within the fixed event windows')
a.set_title('a  GRAFT decoded alternatives',loc='left',weight='bold',pad=13)
a.legend(loc='center left',bbox_to_anchor=(0,.28),frameon=False,fontsize=8,ncol=1)
a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
a=axs[0,1];cs=[comp['comparisons'][k] for k in ['slap_sweep','native_slap']]
for i,d in enumerate(cs):
 n=d['union_classes'];total=d['total_classes'];a.barh(i,n,color=GREEN,height=.43);a.barh(i,total-n,left=n,color=RED,height=.43)
 a.text(n/2,i,f'{n}/{total}',ha='center',va='center',color='white',fontsize=10,weight='bold')
 a.text(0,i-.34,['Against SLAP sweep','Against native SLAP'][i],fontsize=8.5,color=INK)
a.set(yticks=[],xlim=(0,180),xticks=[0,50,100,150],xlabel='Comparator minimum-event classes recovered')
a.set_title('b  GRAFT coverage of SLAP alternatives',loc='left',weight='bold',pad=13);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
a.text(.02,.03,'All compared classes covered in 139/140 reactions.\nRed: classes absent from GRAFT.',transform=a.transAxes,fontsize=8.2,color=MUTED,va='bottom');a.set_ylim(1.85,-.6)
a=axs[1,0];old=np.array([r['old_literal_branches'] for r in rows]);new=np.array([r['final_branches'] for r in rows]);a.plot([10,8e4],[10,8e4],color=LIGHT,ls='--',lw=1)
a.scatter(old,new,s=20,c=GREEN,alpha=.65,edgecolors='white',linewidths=.4)
r=next(r for r in rows if r['case']==25);a.scatter([r['old_literal_branches']],[r['final_branches']],s=42,c=ORANGE,zorder=5);a.annotate('52,669 → 2,856',xy=(52669,2856),xytext=(1000,100),arrowprops=dict(arrowstyle='-',color=ORANGE),color=ORANGE,fontsize=8.5)
a.set(xscale='log',yscale='log',xlim=(12,8e4),ylim=(12,8e4),xlabel='Ordered GRAFT branches / reaction',ylabel='Unordered GRAFT branches / reaction');a.set_title('c  GRAFT fragment deduplication',loc='left',weight='bold',pad=13);a.grid(alpha=.12)
a.text(.03,.97,f'Total: {dedup["old_branches"]:,} → {dedup["new_branches"]:,}\nMedian: {dedup["old_median"]:g} → {dedup["new_median"]:g}',transform=a.transAxes,va='top',fontsize=8.5)
a=axs[1,1];times=[dedup['wall_seconds']/60,dedup['cpu_seconds']/60]
a.barh([0,1],times,color=[GREEN,BLUE],height=.42)
for i,v in enumerate(times):
 a.text(v+.4,i,f'{v:.2f} min',va='center',fontsize=9,weight='bold')
 a.text(0,i-.34,['Wall time (shared load)','CPU time (sum over workers)'][i],fontsize=8.5,color=INK)
a.set(yticks=[],xlim=(0,34),ylim=(1.85,-.6),xticks=[0,10,20,30],xlabel='Complete catalogue and event-window decoding')
a.set_title('d  GRAFT postprocessing cost',loc='left',weight='bold',pad=13);a.grid(axis='x',alpha=.15);a.set_axisbelow(True)
a.text(.02,.03,f'{dedup["workers"]} workers · {dedup["peak_mib"]/1024:.2f} GiB peak combined memory\nSearch and competition excluded.',transform=a.transAxes,fontsize=8.2,color=MUTED,va='bottom')
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

competitors=read('competitors.json');assert competitors['rescored_with_current_evaluator'] and competitors['unchanged_case_outcomes']
labels={'rxnmapper':'RXNMapper','localmapper':'LocalMapper','chython':'Chython','slap_binary':'SLAP, binary','slap_weighted':'SLAP, weighted','indigo':'Indigo','rdt':'RDT'}
rows=[]
def comparison_row(name,setting,output,n):
 return f"{name} & {setting} & {output} & {n:,} ({100*n/N:.2f}\\%) " + r"\\"
for k,nseed in zip(keys,[1,2,3,10]):
 rows.append(comparison_row('GRAFT',f'{nseed} seed'+('s' if nseed>1 else '')+', sweep'+(' (default)' if nseed==1 else ''),'Families',methods[k]['golden_outcomes']['recovered']))
rows.append(comparison_row('GRAFT','1 seed, no sweep (ablation)','Families',uncut['methods']['graft']['counts']['recovered']))
rows.append(r"\midrule")
for d in competitors['methods']:
 assert d['total']==N and len(d['any_correct_cases'])==d['any_correct']
 if not d['method'].startswith('slap'):continue
 rows.append(comparison_row(labels[d['method']],'Default, no sweep','Candidates',d['any_correct']))
rows.append(comparison_row('SLAP, union','No sweep','Candidates',uncut['methods']['slap']['counts']['recovered']))
rows.append(comparison_row('SLAP, union','Sweep','Candidates',slap['sweep_union_recovered']))
rows.append(r"\midrule")
for d in competitors['methods']:
 if d['method'].startswith('slap'):continue
 rows.append(comparison_row(labels[d['method']],'Default, no sweep','Mapping',d['any_correct']))
(MAN/'includes/generated-competitor-table.tex').write_text(r"\begin{tabular}{@{}lllr@{}}\toprule"+'\n'+r"Method & Search setting & Evaluated output & Reference recovered\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')

# The case-level failure table is generated from the audited final misses.
misses=read('golden_miss_analysis.json')
labels={'extra_matched_pair':'One extra matched pair',
 'atom_correspondence_excluded_by_relaxed_orbits':'Atom pairing excluded by orbit bounds',
 'joint_correspondence_excluded_by_full_verifier':'Joint C/O correspondence excluded'}
rows=[f"{r['case']} & {r['reference_pairs']} & {r['retained_pair_counts'][0]} & {labels[r['classification']]} " + r"\\" for r in misses['rows']]
(MAN/'includes/generated-miss-table.tex').write_text(r"\begin{tabular}{lrrl}\toprule"+'\n'+r"Case & Reference pairs & Retained pairs & Observed mismatch\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')

# Direction policies use explicit-atom counts and a stated tie rule.
direction=read('direction_recovery.json');rows=[]
for key,label in [('seeds1','1 seed + sweep (default)'),('seeds2','2 seeds + sweep')]:
 d=direction['methods'][key]['groups']['all'];assert d['n']==N
 cells=[f"{d[k]['recovered']:,} ({d[k]['percent']:.2f}\\%)" for k in ['smaller_first','larger_first','bidirectional']]
 rows.append(label+' & '+' & '.join(cells)+r"\\")
(MAN/'includes/generated-direction-table.tex').write_text(r"\begin{tabular}{@{}lrrr@{}}\toprule"+'\n'+r"GRAFT setting & Smaller-first & Larger-first & Bidirectional\\\midrule"+'\n'+'\n'.join(rows)+'\n'+r"\bottomrule\end{tabular}"+'\n')
