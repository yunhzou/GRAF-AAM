#!/usr/bin/env python3
"""Build the README coverage figure/table from the manuscript's audited evidence."""
import hashlib,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


def build(root):
    output=root/'docs/assets';output.mkdir(parents=True,exist_ok=True)
    evidence=root/'manuscript/evidence'
    names=['seed_comparison.json','slap_sweep.json','competitors.json','unswept.json']
    records={name:json.loads((evidence/name).read_text()) for name in names}
    seeds=records[names[0]]['methods'];slap=records[names[1]];competitors=records[names[2]]
    n=competitors['denominator'];assert n==1851 and competitors['rescored_with_current_evaluator']
    rows=[]
    def row(method,setting,kind,count,key):
        assert 0<=count<=n
        rows.append(dict(method=method,setting=setting,output=kind,recovered=count,denominator=n,coverage_percent=100*count/n,key=key))
    for seed in [1,3,10]:
        data=seeds[f'seeds{seed}'];count=data['golden_outcomes']['recovered']
        assert data['golden_cases']==n and sum(x['outcome']=='recovered' for x in data['per_case'])==count
        row('GRAFT',f'{seed} seed'+('s' if seed>1 else '')+', sweep'+(' (default)' if seed==1 else ''),'Compressed families',count,f'graft{seed}')
    assert slap['cases']==n and sum(x['outcome']=='recovered' for x in slap['per_case'])==slap['sweep_union_recovered']
    row('SLAP','Bidirectional + our sweep','Multiple candidates',slap['sweep_union_recovered'],'slap_sweep')
    row('SLAP','Bidirectional, no sweep','Multiple candidates',records['unswept.json']['methods']['slap']['counts']['recovered'],'slap_bidirectional')
    methods={m['method']:m for m in competitors['methods']}
    for key in ['slap_binary','localmapper','rxnmapper','chython','rdt','indigo']:
        m=methods[key];assert m['total']==n and len(m['any_correct_cases'])==m['any_correct']
        row('SLAP' if key=='slap_binary' else m['name'],'Default, no sweep'+(' (binary)' if key=='slap_binary' else ''),'Multiple candidates' if key=='slap_binary' else 'One bijection',m['any_correct'],key)
    summary=dict(denominator=n,metric='Reference inclusion among returned alternatives; single-bijection accuracy for single-output methods.',graft_branch_cap=100,graft_directions='bidirectional',rows=rows,sources=[dict(path='manuscript/evidence/'+name,sha256=hashlib.sha256((evidence/name).read_bytes()).hexdigest()) for name in names])
    (output/'golden-coverage.json').write_text(json.dumps(summary,indent=2)+'\n')
    table='| Method | Search setting | Output | References covered | Coverage |\n|---|---|---|---:|---:|\n'
    for r in rows:
        bold=r['key']=='graft1';a='**' if bold else ''
        table+=f"| {a}{r['method']}{a} | {a}{r['setting']}{a} | {r['output']} | {a}{r['recovered']:,} / {n:,}{a} | {a}{r['coverage_percent']:.2f}%{a} |\n"
    readme=root/'README.md';text=readme.read_text();start='<!-- golden-coverage-table:start -->';end='<!-- golden-coverage-table:end -->'
    assert start in text and end in text
    readme.write_text(text.split(start)[0]+start+'\n\n'+table+'\n'+end+text.split(end)[1])
    selected=[next(r for r in rows if r['key']==k) for k in ['graft1','graft10','slap_sweep','localmapper','rxnmapper']]
    colors=['#17745b','#58a487','#aaa18a','#869daf','#bac6ce'];labels=['GRAFT · 1 seed (default)','GRAFT · 10 seeds','SLAP + our sweep','LocalMapper','RXNMapper']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(12,6.5),facecolor='#fafbf8');ax=fig.add_axes([.285,.17,.53,.52]);ax.set_facecolor('#fafbf8')
    fig.text(.055,.92,'Golden benchmark',fontsize=24,fontweight='bold',color='#193d33')
    fig.text(.055,.867,'Reference coverage across all 1,851 reactions',fontsize=13,color='#697d71')
    fig.text(.64,.915,'99.08%',fontsize=31,fontweight='bold',color='#17745b');fig.text(.642,.867,'GRAFT · default',fontsize=11,color='#697d71')
    fig.text(.83,.915,'99.41%',fontsize=31,fontweight='bold',color='#58a487');fig.text(.833,.867,'GRAFT · 10 seeds',fontsize=11,color='#697d71')
    for i,(r,c) in enumerate(zip(selected,colors)):
        v=r['coverage_percent'];ax.barh(i,100,height=.49,color='#e8ede6',zorder=1);ax.barh(i,v,height=.49,color=c,zorder=2)
        ax.text(103,i+.015,f'{v:.2f}%',va='center',fontsize=13,fontweight='bold',color='#254c3e')
        ax.text(103,i+.25,f"{r['recovered']:,} / {n:,}",va='center',fontsize=9,color='#7c8e80')
    ax.set_yticks(range(5),labels);ax.tick_params(axis='y',length=0,pad=16,labelsize=11);ax.tick_params(axis='x',length=0,colors='#809184',labelsize=10,pad=9)
    ax.set_xlim(0,100);ax.set_ylim(4.6,-.6);ax.set_xticks([0,25,50,75,100]);ax.xaxis.set_major_formatter(PercentFormatter());ax.grid(axis='x',color='#dce3d8',linewidth=.7,zorder=0);ax.set_axisbelow(True)
    for s in ax.spines.values():s.set_visible(False)
    fig.text(.055,.068,'GRAFT / SLAP: coverage among alternatives. LocalMapper / RXNMapper: one-bijection accuracy.',fontsize=10,color='#697d71')
    fig.text(.055,.032,'Primary strict re-evaluation · GRAFT: bidirectional sweep, cap 100 · all failures and unresolved cases included',fontsize=9,color='#849282')
    for ext in ['png','svg']:fig.savefig(output/f'golden-coverage.{ext}',dpi=180,facecolor=fig.get_facecolor())
    svg=output/'golden-coverage.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig);print('Verified',len(rows),'table rows against archived evidence; rendered Golden coverage figure.')

if __name__=='__main__':build(Path(__file__).resolve().parents[1])
