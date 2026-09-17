import os
import json,gzip,sys,copy,importlib.util,concurrent.futures,time
from pathlib import Path
import numpy as np
REPO=Path(os.environ.get('GRAFT_AUDIT_REPO', Path(__file__).resolve().parents[3]))
sys.path.insert(0,str(REPO/'src'))
from graft.viewers import comparison_document,reaction_html
from graft.alignment.interpolation import internal_coordinate_interpolation
ROOT=Path(os.environ.get('GRAFT_AUDIT_WORKSPACE', Path.cwd()));AUDIT=ROOT/'outputs/chirality-strict-audit-20260917';OUT=AUDIT/'viewer-build'
def dump(p,data):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(data,separators=(',',':')))
def one(i):
    r=json.load(open(AUDIT/f'case{i}.json'));inp=json.load(open(ROOT/f'work/full140_inputs/{i}/input.json'))
    old=json.load(open(ROOT/f'work/chirality-comparison/old-{r["name"]}.json'))
    prior=json.load(gzip.open(ROOT/f'outputs/chirality-holdout-minimum-20260917/case{i}/viewer-data.json.gz','rt'))
    saved=json.load(open(ROOT/f'outputs/final-end-to-end-optimized-20260916/case{i}/candidates.json'));cs={c['id']:c for c in saved['candidates']}
    n=len(inp['reactant']['elements']);rawR=np.array(inp['reactant']['coordinates']);rawP=np.array(inp['product']['coordinates']);rw=np.array(inp['reactant']['wbo']);pw=np.array(inp['product']['wbo']);rb={(a,b) for a in range(n) for b in range(a+1,n) if rw[a,b]>=.2}
    mechs=[];stats=[]
    for m,o in zip(old['mechanisms'],r['old']):
        m=copy.deepcopy(m);m['label']=f"Old {m['id']} · {len(o['geometry']['violations'])} full-orientation reversals";m['strict_audit']=o['geometry'];mechs.append(m)
    for row in r['rows']:
        c=cs[row['id']];result=row['result'];allowed=result['status']=='allowed';mapping={int(a):b for a,b in (result['mapping'] if allowed else c['mapping']).items()};label=f"C{row['candidate']+1} · "+('STRICT VERIFIED' if allowed else 'NO STRICT SOLUTION · raw display only')
        events=[dict(kind=k,r=[a,b],p=[mapping[a],mapping[b]],wbo=[float(rw[a,b]),float(pw[mapping[a],mapping[b]])]) for k,bs in c['events'].items() for a,b in bs]
        doc=comparison_document(dict(index=i,name=r['name'],endpoints=[inp['reactant'],inp['product']],records=[dict(name=label,mapping=sorted(mapping.items()),events=events,pattern=c['id'])]))
        m=doc['mechanisms'][0]
        cached=next((p for p in [*old['mechanisms'],*prior['mechanisms']] if {int(a):b for a,b in p['mapping_RP'].items()}==mapping),None)
        if cached:path=cached['endpoint_interpolation']
        else:
            pb={(a,b) for a in range(n) for b in range(a+1,n) if pw[mapping[a],mapping[b]]>=.2}
            path=internal_coordinate_interpolation(rawR,rawP[[mapping[a] for a in range(n)]],inp['reactant']['elements'],bonded_pairs=rb|pb,persistent_bonded_pairs=rb&pb,reactant_bonded_pairs=rb,product_bonded_pairs=pb,n_frames=101)
        m['endpoint_interpolation']=path;m['index_chirality']=result['diagnostics'];m['strict_status']=result['status'];mechs.append(m)
        counts=[f['clashes']['count'] for f in path['frames']];stats.append(dict(candidate=row['candidate'],status=result['status'],peak_clashes=max(counts),sum_frame_clashes=sum(counts),clashing_frames=sum(c>0 for c in counts)))
    good=sum(x['result']['status']=='allowed' for x in r['rows']);total=len(r['rows']);tag='ALL STRICT' if good==total else 'NO STRICT SOLUTION' if good==0 else 'MIXED'
    step=f'{i:03d} · {r["name"]} [{good}/{total} strict]'
    note=f'{good}/{total} current minimum-event candidates have a strict orientation-preserving mapping. '
    note+='NO STRICT SOLUTION means none exists in the current saved AAM families with that candidate’s concrete bond events; its displayed raw mapping is NOT a corrected solution. '
    note+='Old mappings are reference displays, with all defined orientation reversals counted independently. Preserving these signs is not a chemical-pathway or clash-free guarantee. '
    note+='101-frame interpolation; magenta marks nonbonded distances < 0.70 × summed covalent radii.'
    doc=comparison_document(dict(index=i,name=step,endpoints=[inp['reactant'],inp['product']],records=[],note=note));doc['mechanisms']=mechs
    for j,m in enumerate(mechs,1):m['id']=j
    doc['default_mech_id']=next((m['id'] for m in mechs if m.get('strict_status')=='allowed'),1)
    dump(OUT/'cases'/step/'summary.json',dict(elapsed_seconds=sum(row['wall_seconds'] for row in r['rows'])))
    dump(OUT/'cases'/step/'rp_stage.json',dict(mechanisms=mechs))
    path=OUT/'views'/step/'view.html';path.parent.mkdir(parents=True,exist_ok=True);path.write_text(reaction_html(doc))
    dump(AUDIT/f'viewer-stats{i}.json',stats)
    return dict(source_index=i,step_id=step,atom_count=n),len(mechs)
if __name__=='__main__':
    specs={k:[] for k in ('small','medium','large')};count=0
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as ex:
        for spec,n in ex.map(one,range(140)):
            specs['small' if spec['atom_count']<50 else 'medium' if spec['atom_count']<100 else 'large'].append(spec);count+=n
    for tier,rows in specs.items():dump(OUT/'manifests'/f'{tier}.json',dict(cases=rows))
    dump(OUT/'inventory.json',dict(case_count=140));dump(OUT/'batch_summary.json',dict(error_count=0,total_elapsed_case_seconds=sum(sum(row['wall_seconds'] for row in json.load(open(AUDIT/f'case{i}.json'))['rows']) for i in range(140))))
    spec=importlib.util.spec_from_file_location('integrated',REPO/'tools/build_integrated_batch_viewer.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    result=mod.build(OUT,AUDIT/'viewer.html',compact=True)
    path=AUDIT/'viewer.html'
    text=path.read_text().replace('AAM integrated viewer — 140 cases','Strict orientation audit — 140 cases')
    text=text.replace(' mechanisms</span>',' mapping views</span>')
    text=text.replace('${c.mechanisms} mechanisms, violations=${c.violations}','${c.mechanisms} mapping views')
    text=text.replace(" mechanism${c.mechanisms===1?'':'s'}", " mapping view${c.mechanisms===1?'':'s'}")
    path.write_text(text)
    print('cases',140,'mapping views',count,'output',path)

